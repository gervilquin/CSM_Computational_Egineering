function [K,Kb,Ka,Ks,Kt,Bb,Ba,Bs,Bt,M,Me,R,Nek,le,w] = ComputeKMmatricesBeam(Ndof,X,Tn,Tm,j_p,E,A,Iy,Iz,G,ky,kz,kt,J,rho)

    % Variables
    nel = size(Tn,1);

    % Initialization
    K = sparse(Ndof,Ndof);
    M = sparse(Ndof,Ndof);
    R = zeros(2*6,2*6,nel);

    Bb = zeros(2,2*6,nel);
    Bs = zeros(2,2*6,nel);
    Ba = zeros(1,2*6,nel);
    Bt = zeros(1,2*6,nel);

    Kb = zeros(2*6,2*6,nel);
    Ks = zeros(2*6,2*6,nel);
    Ka = zeros(2*6,2*6,nel);
    Kt = zeros(2*6,2*6,nel);
    Nek = zeros(6,12,nel,2);
    Me = zeros(2*6,2*6,nel);
    le = zeros(nel,1);

    % Assembly process
    for e = 1:nel
        % Compute rotation matrix
        le(e) = norm(X(Tn(e,2),:) - X(Tn(e,1),:));
        i_p = (X(Tn(e,2),:)-X(Tn(e,1),:))/le(e);
        j_p = j_p(Tm(e),:); 
        k_p = cross(i_p,j_p);

        R_prime = [i_p' j_p' k_p'      zeros(3);
                     zeros(3)     i_p' j_p' k_p']';
        R(:,:,e) = [R_prime   zeros(6);
                    zeros(6) R_prime  ];
        
        % Compute the shape function derivates
        Nx = [-1/le(e) 1/le(e)];

        %Compute each element matrix
            % Axial component of stiffness matrix
            Ba(1,:,e) = [Nx(1) 0 0 0 0 0 Nx(2) 0 0 0 0 0];
            Ca = E(Tm(e))*A(Tm(e)); % Only one material
            Ka(:,:,e) = le(e)*R(:,:,e)'*Ba(:,:,e)'*Ca*Ba(:,:,e)*R(:,:,e);

            %Bending component of stiffness matrix
            Bb(:,:,e) = [0 0 0 0 Nx(1)   0   0 0 0 0 Nx(2)   0  ;
                         0 0 0 0   0   Nx(1) 0 0 0 0   0   Nx(2)] ;
            Cb = E(Tm(e))*[Iy(Tm(e))     0    ;
                               0     Iz(Tm(e))];
            Kb(:,:,e) = le(e)*R(:,:,e)'*Bb(:,:,e)'*Cb*Bb(:,:,e)*R(:,:,e);

            % Shear component of stiffness matrix
            Ns = 1/2;
            Bs(:,:,e) = [0 Nx(1)  0    0 0 -Ns 0 Nx(2)  0    0 0 -Ns;
                         0  0    Nx(1) 0 Ns 0  0  0    Nx(2) 0 Ns 0 ]; 

            Cs = G(Tm(e))*A(Tm(e))*[ky(Tm(e))     0    ; 
                                        0     kz(Tm(e))];
            Ks(:,:,e) = le(e)*R(:,:,e)'*Bs(:,:,e)'*Cs*Bs(:,:,e)*R(:,:,e);

            % Torsion component of stiffness matrix
            Bt(:,:,e) = [0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0 0];
            
            Ct = G(Tm(e))*J(Tm(e))*kt(Tm(e));

            Kt(:,:,e) = le(e)*R(:,:,e)'*Bt(:,:,e)'*Ct*Bt(:,:,e)*R(:,:,e);

            % Mass matrix
            xi = [-1/sqrt(3), 1/sqrt(3)]; % gauss points coordinates
            w = [1 1]; % gauss points weights

            rho_matrix = zeros(6);
            rho_matrix(1,1) = A(Tm(e)); rho_matrix(2,2) = A(Tm(e)); rho_matrix(3,3) = A(Tm(e)); rho_matrix(4,4) = J(Tm(e)); rho_matrix(5,5) = Iy(Tm(e)); rho_matrix(6,6) = Iz(Tm(e));
            rho_matrix = rho(Tm(e))*rho_matrix;

            for k=1:2
                N_1 = (1-xi(k))/2;
                N_2 = (1+xi(k))/2;
                Nek(:,:,e,k) = [N_1*eye(6) N_2*eye(6)];
                Me(:,:,e) = Me(:,:,e) + w(k)*le(e)*R(:,:,e)'*Nek(:,:,e,k)'*rho_matrix*Nek(:,:,e,k)*R(:,:,e)/2;
            end

        % Assembly to global matrices
        Idof = zeros(12,1);
        for j =1:6
            Idof(j) = 6*(Tn(e,1)-1) + j;
            Idof(6+j) = 6*(Tn(e,2)-1) + j;
        end

        K(Idof,Idof) = K(Idof,Idof) + Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e);
        M(Idof,Idof) = M(Idof,Idof) + Me(:,:,e);
    
    end
end