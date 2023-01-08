function [F] = ComputeFvector(Ne,X,Tn,Pe)
    % Initialization
    nDOFs = ((Ne+1)^2)*6;
    F = zeros(nDOFs,1);

    % Nodal distributed forces vector
    P = zeros((Ne+1)^2,6);
    for i=1:length(Pe)
        P(Pe(i,2),Pe(i,3)) = Pe(i,1);
    end
    
    p = zeros(4*6,1);
    N_mass = zeros(5,4*5,Ne^2,4);
    Fe = zeros(4*6,Ne^2);
    for e=1:Ne^2
        p(:,e) = [P(Tn(e,1),:) P(Tn(e,2),:) P(Tn(e,3),:) P(Tn(e,4),:)]';

        % Rotation matrix
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

        % Shape functions computation (nG = 4)
        nG = 4;
        [N,derN,wk] = ComputeShapeFunction(nG);
        N_node = zeros(5,5,4);
        for i=1:4
            N_node(:,:,i) = N(i)*eye(5);
        end
        
        for j=1:nG
            J = zeros(2);  
            for i=1:4 % 4 nodes for quadrilateral elements
                J = J + [derN(i,1,j) derN(i,2,j)]'*X(Tn(e,i),:)*[i_prime j_prime];
            end
            S(e,j) = wk(j)*det(J);
            N_mass(:,:,e,j) = [N_node(:,:,1) N_node(:,:,2) N_node(:,:,3) N_node(:,:,4)];

            % Element force vector computation
            Fe(:,e) = Fe(:,e) + S(e,j)*R(:,:,e)'*N_mass(:,:,e,j)'*N_mass(:,:,e,j)*R(:,:,e)*p(:,e);
        end
        
        % Global force vector assembly
        Idof = zeros(24,1);
        for j=1:6
            Idof(j,1) = 6*(Tn(e,1)-1)+j;
            Idof(6+j,1) = 6*(Tn(e,2)-1)+j;
            Idof(12+j,1) = 6*(Tn(e,3)-1)+j;
            Idof(18+j,1) = 6*(Tn(e,4)-1)+j;       
        end
        F(Idof,1) = F(Idof,1) + Fe(:,e);
    end
end