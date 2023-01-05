function [KB,KBb,KBa,KBs,KBt,BBb,BBa,BBs,BBt,MB,MBe,RB,NekB,leB,xi,w] = beam_global_matrices_assembly(Nnod,Ndof,X,Tn,Tm,j_p,E,A,Iy,Iz,G,ky,kz,kt,J,rho)

    % Initialization
    KB = sparse(Ndof,Ndof);
    KBb = zeros(12,12,Nnod);
    KBs = zeros(12,12,Nnod);
    KBa = zeros(12,12,Nnod);
    KBt = zeros(12,12,Nnod);

    BBb = zeros(2,12,Nnod);
    BBs = zeros(2,12,Nnod);
    BBa = zeros(1,12,Nnod);
    BBt = zeros(1,12,Nnod);

    MB = sparse(Ndof,Ndof);
    MBe = zeros(12,12,Nnod);
    RB = zeros(12,12,Nnod);
    NekB = zeros(6,12,Nnod,2);
    leB = zeros(Nnod,1);

    % Assembly process
    for e = 1:length(Tn)
        % Compute rotation matrix
        leB(e) = norm(X(Tn(e,2),:) -X(Tn(e,1),:));
        i_p = (X(Tn(e,2),:)-X(Tn(e,1),:))/leB(e);
        j_p = j_p(Tm(e),:); 
        k_p = cross(i_p,j_p);

        R_node = [i_p' j_p' k_p' zeros(3);
            zeros(3) i_p' j_p' k_p']';
        RB(:,:,e) = [R_node zeros(6);
            zeros(6) R_node];
        
        % Compute the shape function derivates
        Nx = [-1/leB(e) 1/leB(e)];

        %Compute each element matrix
            % Axial component of stiffness matrix
            BBa(:,:,e) = [Nx(1) 0 0 0 0 0 Nx(2) 0 0 0 0 0];
            Ca = E(Tm(e))*A(Tm(e)); % Only one material
            KBa(:,:,e) = leB(e)*RB(:,:,e)'*BBa(:,:,e)'*Ca*BBa(:,:,e)*RB(:,:,e);

            %Bending component of stiffness matrix
            BBb(:,:,e) = [0 0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0;
                    0 0 0 0 0 Nx(1) 0 0 0 0 0 Nx(2)] ;
            Cb = E(Tm(e))*[Iy(Tm(e)) 0;
                        0 Iz(Tm(e))];
            KBb(:,:,e) = leB(e)*RB(:,:,e)'*BBb(:,:,e)'*Cb*BBb(:,:,e)*RB(:,:,e);

            % Shear component of stiffness matrix
            Ns = 1/2;
            BBs(:,:,e) = [0 Nx(1) 0 0 0 -Ns 0 Nx(2) 0 0 0 -Ns;
                  0 0 Nx(1) 0 Ns 0 0 0 Nx(2) 0 Ns 0]; 

            Cs = G(Tm(e))*A(Tm(e))*[ky(Tm(e)) 0; 
                        0 kz(Tm(e))];
            KBs(:,:,e) = leB(e)*RB(:,:,e)'*BBs(:,:,e)'*Cs*BBs(:,:,e)*RB(:,:,e);

            % Torsion component of stiffness matrix
            BBt(:,:,e) = [0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0 0];
            
            Ct = G(Tm(e))*J(Tm(e))*kt(Tm(e));

            KBt(:,:,e) = leB(e)*RB(:,:,e)'*BBt(:,:,e)'*Ct*BBt(:,:,e)*RB(:,:,e);

            % Mass matrix
            xi = [-1/sqrt(3), 1/sqrt(3)]; % gauss points coordinates
            w = [1 1]; % gauss points weights

            rho_ = zeros(6);
            rho_(1,1) = A(Tm(e)); rho_(2,2) = A(Tm(e)); rho_(3,3) = A(Tm(e)); rho_(4,4) = J(Tm(e)); rho_(5,5) = Iy(Tm(e)); rho_(6,6) = Iz(Tm(e));
            rho_ = rho(Tm(e))*rho_;

            %Me = zeros(12);
            
            for k=1:2
                N_1 = (1-xi(k))/2;
                N_2 = (1+xi(k))/2;
                NekB(:,:,e,k) = [N_1*eye(6) N_2*eye(6)];
                MBe(:,:,e) = MBe(:,:,e) + w(k)*leB(e)*RB(:,:,e)'*NekB(:,:,e,k)'*rho_*NekB(:,:,e,k)*RB(:,:,e)/2;
            end

        % Assembly to global matrices
        Idof = zeros(1,12);
        for j =1:6
            Idof(j) = 6*(Tn(e,1)-1) + j;
            Idof(6+j) = 6*(Tn(e,2)-1) + j;
        end

        KB(Idof,Idof) = KB(Idof,Idof) + KBa(:,:,e) + KBb(:,:,e) + KBs(:,:,e) + KBt(:,:,e);
        MB(Idof,Idof) = MB(Idof,Idof) + MBe(:,:,e);
    
    end % end loop element
    



end