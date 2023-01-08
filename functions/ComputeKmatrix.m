function [K, Bs, Bmt, Bmn, Bb, R] = ComputeKmatrix(Ne,X,Tn,E,h,nu)

    % Variables
    nDOFs = ((Ne+1)^2)*6;
    
    % Initialisation
    K = zeros(nDOFs);
%     M = zeros(nDOFs);
    R = zeros(4*5,4*6,Ne^2);

    Bs = zeros(2,4*5,Ne^2);
    Bmt = zeros(1,4*5,Ne^2);
    Bmn = zeros(2,4*5,Ne^2,4);
    Bb = zeros(3,4*5,Ne^2,4);

    Ks = zeros(4*6,4*6,Ne^2);
    Km = zeros(4*6,4*6,Ne^2);
    Kb = zeros(4*6,4*6,Ne^2);
%     N_mass = zeros(5,4*5,Ne^2,4);
%     Me = zeros(4*6,4*6,Ne^2);

    % Assembly process
    for e=1:Ne^2
        % rotation matrix
        S = cross((X(Tn(e,3),:)'-X(Tn(e,1),:)'),(X(Tn(e,4),:)'-X(Tn(e,2),:)')/2);
        k_prime = S/norm(S);
        d = (X(Tn(e,2),:)'-X(Tn(e,1),:)'+X(Tn(e,3),:)'-X(Tn(e,4),:)')/2;
        i_prime = d/norm(d);
        j_prime = cross(k_prime,i_prime);
        R_prime = [i_prime j_prime k_prime      zeros(3,2);
                      zeros(3,3)            i_prime j_prime]';
        R(:,:,e) = [  R_prime    zeros(5,6)   zeros(5,6)   zeros(5,6);
                     zeros(5,6)   R_prime     zeros(5,6)   zeros(5,6);
                     zeros(5,6)  zeros(5,6)    R_prime     zeros(5,6);
                     zeros(5,6)   zeros(5,6)   zeros(5,6)   R_prime  ];

        % Shape functions computation (nG = 1)
        nG = 1;
        [N,derN,~] = ComputeShapeFunction(nG); 
        for j=1:nG
            J = zeros(2);  
            for i=1:4 % 4 nodes for quadrilateral elements
                J = J + [derN(i,1,j) derN(i,2,j)]'*X(Tn(e,i),:)*[i_prime j_prime];
            end
            N_global = J\[derN(:,1,j) derN(:,2,j)]';
            S = 4*det(J);

            % Shear stiffness matrix
            Bs_node = zeros(2,5,4);
            for i=1:4
                Bs_node(:,:,i) = [0  0   N_global(1,i)   0     N(i);
                                 0  0   N_global(2,i)  -N(i)    0 ];
            end
            Cs = eye(2)*5*h*E/(12*(1+nu));
            Bs(:,:,e) = [Bs_node(:,:,1) Bs_node(:,:,2) Bs_node(:,:,3) Bs_node(:,:,4)];
            Ks(:,:,e) = S*R(:,:,e)'*Bs(:,:,e)'*Cs*Bs(:,:,e)*R(:,:,e);

            % Membrane transverse stiffness matrix
            Bmt_node = zeros(1,5,4);
            for i=1:4
                Bmt_node(:,:,i) = [N_global(2,i) N_global(1,i)  0   0   0];
            end
            Cmt = h*E/(2*(1+nu));
            Bmt(:,:,e) = [Bmt_node(:,:,1) Bmt_node(:,:,2) Bmt_node(:,:,3) Bmt_node(:,:,4)];
            Km(:,:,e) = S*R(:,:,e)'*Bmt(:,:,e)'*Cmt*Bmt(:,:,e)*R(:,:,e);
        end

        % Shape functions computation (nG = 4)
        nG = 4;
        [~,derN,wk] = ComputeShapeFunction(nG); 
        for j=1:nG
            J = zeros(2);  
            for i=1:4 % 4 nodes for quadrilateral elements
                J = J + [derN(i,1,j) derN(i,2,j)]'*X(Tn(e,i),:)*[i_prime j_prime];
            end
            N_global = J\[derN(:,1,j) derN(:,2,j)]';
            S(e,j) = wk(j)*det(J);
            
            % Membrane normal stiffness matrix
            Bmn_node = zeros(2,5,4);
            for i=1:4
                Bmn_node(:,:,i) = [N_global(1,i)       0        0    0   0;
                                    0             N_global(2,i)  0    0   0];
            end
            Cmn = h*E/(1-nu^2)*[ 1    nu;
                                nu     1];
            Bmn(:,:,e,j) = [Bmn_node(:,:,1) Bmn_node(:,:,2) Bmn_node(:,:,3) Bmn_node(:,:,4)];
            Km(:,:,e) = Km(:,:,e) + S(e,j)*R(:,:,e)'*Bmn(:,:,e,j)'*Cmn*Bmn(:,:,e,j)*R(:,:,e);

            % Bending moment stiffness matrix
            Bb_node = zeros(3,5,4);
            for i=1:4
                Bb_node(:,:,i) =  [0  0   0        0         N_global(1,i);
                                   0  0   0   N_global(2,i)       0       ;
                                   0  0   0  -N_global(1,i)  N_global(2,i)];
                             
            end
            Cb = h^(3)*E/(12*(1-nu^2))*[ 1   nu      0   ;
                                        nu    1      0   ;
                                         0    0  (1-nu)/2];
            Bb(:,:,e,j) = [Bb_node(:,:,1) Bb_node(:,:,2) Bb_node(:,:,3) Bb_node(:,:,4)];
            Kb(:,:,e) = Kb(:,:,e) + S(e,j)*R(:,:,e)'*Bb(:,:,e,j)'*Cb*Bb(:,:,e,j)*R(:,:,e);

%             % Mass matrix
%             N_node = zeros(5,5,4);
%             for i=1:4
%                 N_node(:,:,i) = N(i)*eye(5);
%             end
%             rho_matrix = rho*h*[1 0 0     0      0  ;
%                                 0 1 0     0      0  ;
%                                 0 0 1     0      0  ;
%                                 0 0 0  h^2/12    0  ;
%                                 0 0 0     0   h^2/12];
%             N_mass(:,:,e,j) = [N_node(:,:,1) N_node(:,:,2) N_node(:,:,3) N_node(:,:,4)];
%             Me(:,:,e) = Me(:,:,e) + S(e,j)*R(:,:,e)'*N_mass(:,:,e,j)'*rho_matrix*N_mass(:,:,e,j)*R(:,:,e);
        end

        % Global K matrix assembly
        Idof = zeros(24,1);
        for j=1:6
            Idof(j,1) = 6*(Tn(e,1)-1)+j;
            Idof(6+j,1) = 6*(Tn(e,2)-1)+j;
            Idof(12+j,1) = 6*(Tn(e,3)-1)+j;
            Idof(18+j,1) = 6*(Tn(e,4)-1)+j;
        end
        K(Idof,Idof) = K(Idof,Idof) + Km(:,:,e) + Kb(:,:,e) + Ks(:,:,e);
%         M(Idof,Idof) = M(Idof,Idof) + Me(:,:,e);
    end

    % Artificial rotation stiffness matrix
    n = zeros(3,(Ne+1)^2);
    k = zeros(3,Ne^2);
    for e=1:Ne^2
        % Compute normal and surface
        S_matrix = cross((X(Tn(e,3),:)'-X(Tn(e,1),:)'),(X(Tn(e,4),:)'-X(Tn(e,2),:)')/2);
        S(e) = sqrt(S_matrix(1)^2 + S_matrix(2)^2 + S_matrix(3)^2);
        k(:,e) = S_matrix/S(e);

        % Assemble to get nodal normal
        for i=1:4
            n(:,Tn(e,i)) = n(:,Tn(e,i)) + k(:,e);
        end
    end
    % Compute artificial rotation matrix
    Kr = zeros(nDOFs);
    for e=1:Ne^2
        for i=1:4
            angle = rad2deg( acos( dot(n(:,Tn(e,i)),k(:,e))/sqrt(dot(n(:,Tn(e,i)),n(:,Tn(e,i)))) ) );
            if angle < 5
                Idof = 6*(Tn(e,i)-1) + [4,5,6]';
                Kr(Idof,Idof) = Kr(Idof,Idof) + E*h*S(e)*k(:,e)*k(:,e)';
            end
        end
    end
    K = K + Kr;
end
