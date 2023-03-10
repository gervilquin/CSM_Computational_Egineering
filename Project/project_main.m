clc
clear
close all
addpath("functions\")

JumpToSolver = true; % Set to true once you are confident about the
                      % precomputations (mass and stiffness matrices and 
                      % force vector assemblies) to go directly to the 
                      % solver parts

if ~JumpToSolver

%% 1) Input data

% Mesh information
load('InputData.mat','X','Tn_s','Tn_b','Tm_s','Tm_b','I_root','n_u','n_l','I_le','I_te');

Nnod = length(X(:,1));
Ndof = Nnod*6;
AoA = deg2rad(10);
p_inf = 1.5625e5;
%p_inf = 0;
g = -9.81;

% Material properties
    %1. Skin
    %2. Spars
    %3. Ribs
    %4. Stringers

    % Youngs modulus
    % material: 1       2       3       
    ES = [       200     180     100]*(10^9); %Pa
    EB = [       190]*(10^9); %Pa

    % Poisson ratio
    nuS = [      0.27    0.3     0.33];
    nuB = [      0.3];

    % Density
    rhoS = [     1500    3300    2000]; %kg/m^3
    rhoB = [     2200]; %kg/m^3

    % Thickness / diameter
    hS = [       1.5     15      3]*10^(-3); %m
    dB = [       10]*10^(-3); %m

% Prescribed degrees of freedom
Up = AssignPrescribedDisplacements(0,X); % Prescribed displacements are 0 for all DOFs of each prescribed node

% Point forces (N)
Fe = [0 1 1];
% Body forces (N/kg) -- weight
Be = ComputeBodyForces(X,g);
% Distributed surface forces (N/m2) -- pressure
PeB = [0 1 1];
PeS = ComputeSurfaceForces(X,p_inf,n_u,n_l,AoA);

plot_pressuredistribution(X,n_l,n_u,p_inf,AoA)
%% BEAMS matrices

% Compute beam section properties
[GB,j_pB,AB,JB,IyB,IzB,kyB,kzB,ktB] = BeamSectionProperties(EB(1),nuB(1),dB(1));
%[G,j_p,A,J,Iy,Iz,ky,kz,kt] = BeamSectionProperties(EB(1),nuB(1),dB(1));
%GB = [0 0 0 G];
% GB = [G];
% j_pB = [j_p];
% AB = [A]; JB = [J];
% IyB = [Iy]; IzB = [Iz];
% kyB = [ky]; kzB = [kz]; ktB = [kt];

% Compute stiffness and mass matrix
[KB,KBb,KBa,KBs,KBt,BBb,BBa,BBs,BBt,MB,MBe,RB,NekB,leB,w] = ComputeKMmatricesBeam(Ndof,X,Tn_b,Tm_b,j_pB,EB,AB,IyB,IzB,GB,kyB,kzB,ktB,JB,rhoB);
% Compute force vector
[FB] = ComputeFvectorBeam(Ndof,Tn_b,Fe,PeB,Be,MBe,RB,w,NekB,leB,true);
                         
%% SHELLS matrices

% Compute stiffness and mass matrix
[KS, MS, BSs, BSmt, BSmn, BSb, MSe, RS, S_4, NS_mass] = ComputeKMmatricesShell(Ndof,X,Tn_s,Tn_b,Tm_s,ES,hS,nuS,rhoS);
%% Compute force vector 
[FS] = ComputeFvectorShell(Ndof,Tn_s,Fe,PeS,Be,MSe,RS,S_4,NS_mass,false);

%% Boundary conditions
[If, Ip, u] = ComputeBoundaryConditions(Ndof,Up);

%% Save data
save('Variables.mat');
else    
% Load precompued data
load('Variables.mat');
end

%% Solve system
% Assembly mass and stifness matrix
K = KB + KS;
M = MB + MS;
F = FB + FS;

% Equations solver
[u,FR] = SystemSolver(K,F,u,If,Ip);

%% Compute strain and stress

% (NOT NECESSARY UNTIL NOW) [Sa,Ss,St,Sb,Fx,Fy,Fz,Mx,My,Mz] = compute_interal_forces_strain(Nnod-1,Tn_b,u,BBa,BBs,BBt,BBb,RB,KBa,KBb,KBs,KBt);
[sigVM] = ComputeVonMissesStresses(Tn_s, Tm_s, u, BSs, BSmt, BSmn, BSb, RS, nuS, ES, hS);


%% Plot (a) - deformed state and stress distribution

scale = 5; % Set appropriate scale to visualize the deformation
plotWing(X,Tn_s,Tm_s,u,scale,sigVM);
    
%% Modal analysis

% define the number of modes that want to be returned
Nm = 12; %first 6 modes

% Obtain the first eigenvectos and eigenvalues
[V,D] = eigs(K(If,If),M(If,If),Nm,'sm');

% Obtain the natural frequencies and the vibration modes
Phi = zeros(Ndof,Nm); 
w2 = zeros(1,Nm);

for k =1:length(V(1,:))
    Phi(If,k) = V(:,k)/sqrt(V(:,k)'*M(If,If)*V(:,k));
    w2(k) = D(k,k);
end

%% Plot (b) - vibration modes

plotModes(X,Tn_s,Phi,w2);

%% Reduced order model

% 2-modes
N_modes = 2;
Im = linspace(1,N_modes,N_modes);
W_k = 0; % static case
F_k = F;
u2 = zeros(Ndof,1);

for k =1:length(F_k(1,:))
    for j = 1:length(Im)
        alpha_jk = Phi(:,Im(j))'*F_k(:,k)/(w2(j)-W_k(k)^2);
        u2(:,k) = u2(:,k)+  Phi(:,Im(j))*alpha_jk;
    end
end

[SigVM2] = ComputeVonMissesStresses(Tn_s, Tm_s, u2, BSs, BSmt, BSmn, BSb, RS, nuS, ES, hS);

plotWing(X,Tn_s,Tm_s,u2,scale,SigVM2);

% 6-modes
N_modes = 6;
Im = linspace(1,N_modes,N_modes);
W_k = 0; % static case
F_k = F;
u6 = zeros(Ndof,1);

for k =1:length(F_k(1,:))
    for j = 1:length(Im)
        alpha_jk = Phi(:,Im(j))'*F_k(:,k)/(w2(j)-W_k(k)^2);
        u6(:,k) = u6(:,k)+  Phi(:,Im(j))*alpha_jk;
    end
end

[SigVM6] = ComputeVonMissesStresses(Tn_s, Tm_s, u2, BSs, BSmt, BSmn, BSb, RS, nuS, ES, hS);

plotWing(X,Tn_s,Tm_s,u6,scale,SigVM6);

%% Plot (c) - leading and trailing edges vertical displacement
plot_wing_LE_TE(X,u,u2,u6,I_le,I_te)