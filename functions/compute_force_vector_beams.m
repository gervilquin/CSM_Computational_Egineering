function [F] = compute_force_vector_beams(Ne,Fe,Qe,Be,Tn,Me,R,Nek,le,xi,w)
    
    % Point loads
    N_nodes = Ne+1;
    Ndof = 6*(Ne+1);
    F = zeros(Ndof,1);

    for q = 1:length(Fe(:,1))
        F(6*(Fe(q,2)-1) + Fe(q,3),1) = F(6*(Fe(q,2)-1) + Fe(q,3),1) + Fe(q,1);
    end

    % Nodal distributed loads
    Q = zeros(N_nodes,6);

    for r=1:length(Q(:,1))
        Q(Qe(r,2),Qe(r,3)) = Qe(r,1);
    end

    % Nodal body forces
    B = zeros(Ne+1,6);
    for s=1:length(Be(:,1))
        B(Be(s,2),Be(s,3)) = Be(s,1);
    end

    % Assembly process
    for e =1:length(Tn(:,1))
        % Compute the eleent force vector
        be = [B(Tn(e,1),:),B(Tn(e,2),:)]';
        qe = [Q(Tn(e,1),:),Q(Tn(e,2),:)]';
        fe = Me(:,:,e)*be;
        for k=1:2
            fe = fe+w(k)*le(e)*R(:,:,e)'*Nek(:,:,e,k)'*Nek(:,:,e,k)*R(:,:,e)*qe/2;
        end
        % Assembly to global force vector
        Idof = zeros(1,12);
        for j =1:6
            Idof(j) = 6*(Tn(e,1)-1) + j;
            Idof(6+j) = 6*(Tn(e,2)-1) + j;
        end
        F(Idof,1) = F(Idof,1) + fe; 

    
    end



end