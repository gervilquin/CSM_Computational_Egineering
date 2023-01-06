clear
close all
addpath("functions\")

JumpToSolver = false; % Set to true once you are confident about the
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
g = -9.81;

% Material properties
    %1. Skin
    %2. Spars
    %3. Ribs
    %4. Stringers

    % Youngs modulus
    % material: 1       2       3       
    ES = [       200     180     100]; %GPa
    EB = [       190]; %GPa

    % Poisson ratio
    nuS = [      0.27    0.3     0.33];
    nuB = [      0.3];

    % Density
    rhoS = [     1500    3300    2000]; %kg/m^3
    rhoB = [     2200]; %kg/m^3

    % Thickness / diameter
    hS = [       1.5     15      3]; %mm
    dB = [       10]; %mm

% Prescribed degrees of freedom
Up = prescribed_dof(0,X);

% Distributed body forces (N/kg) -- weight
Be = body_forces(X,g);

% Punctual forces
Fe = [0 1 1];

% Distributed surface forces (N/m2) -- pressure
Pe = surface_forces(X,p_inf,n_u,n_l,AoA);
Qe = [  zeros(Nnod,1),   linspace(1,Nnod,Nnod)',   3*ones(Nnod,1)];

%% BEAMS matrices

% Compute beam section properties
[G,j_p,A,J,Iy,Iz,ky,kz,kt] = beam_section_properties(EB(1),nuB(1),dB(1));
GB = [0 0 0 G];
j_pB = [j_p];
AB = [A];
JB = [J];
IyB = [Iy];
IzB = [Iz];
kyB = [ky];
kzB = [kz];
ktB = [kt];

% Element matrices and assembly
[KB,KBb,KBa,KBs,KBt,BBb,BBa,BBs,BBt,MB,MBe,RB,NekB,leB,xi,w] = beam_global_matrices_assembly(Nnod,Ndof,X,Tn_b,Tm_b,j_p,EB,AB ,IyB,IzB,GB,kyB,kzB,ktB,JB,rhoB);

%% SHELLS matrices

% Initialization
KS = sparse(Ndof,Ndof);
MS = sparse(Ndof,Ndof);

% Element matrices and assembly
% ...

%% Forces vector assembly

f = zeros(Ndof,1);
P = zeros(Nnod,6);
B = zeros(Nnod,6);

% Force vector assembly
[FB] = compute_force_vector_beams(Nnod-1,Fe,Qe,Be,Tn_b,MBe,RB,NekB,leB,xi,w);

%% Boundary conditions

[u,If,Ip] = Compute_boundary_conditions(Nnod-1,Up);

%% Save data

save('Variables.mat');

else
%% Load precompued data

load('Variables.mat');

end

%% Solve system
K = KB;
F = FB;

[u,FR] = solve_system(u,K,F,If,Ip);

%% Compute strain and stress

[Sa,Ss,St,Sb,Fx,Fy,Fz,Mx,My,Mz] = compute_interal_forces_strain(Nnod-1,Tn_b,u,BBa,BBs,BBt,BBb,RB,KBa,KBb,KBs,KBt);

%% Plot (a) - deformed state and stress distribution

scale = 1; % Set appropriate scale to visualize the deformation

% only for testing
SigVM = Sa;
plotWing(X,Tn_s,Tm_s,u,scale,SigVM);
    
%% Modal analysis

% ...

%% Plot (b) - vibration modes

plotModes(X,Tn_s,Phi,w2);

%% Reduced order model

% 2-modes
% ...

plotWing(X,Tn_s,Tm_s,u2,scale,SigVM2);

% 6-modes
% ...

plotWing(X,Tn_s,Tm_s,u6,scale,SigVM6);

%% Plot (c) - leading and trailing edges vertical displacement
% ...