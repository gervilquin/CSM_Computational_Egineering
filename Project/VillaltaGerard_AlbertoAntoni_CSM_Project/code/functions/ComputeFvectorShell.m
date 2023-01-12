function [F] = ComputeFvectorShell(nDOFs,Tn,Fe,PeS,Be,Me,R,S_4,N_mass,FeNotComputed)
    % Initialization
    nel = size(Tn,1);
    nnodes = 4*nel;
    F = zeros(nDOFs,1);

    % Point loads 
    if FeNotComputed == true 
        for i= 1:length(Fe(:,1))
            F(6*(Fe(i,2)-1) + Fe(i,3),1) = F(6*(Fe(i,2)-1) + Fe(i,3),1) + Fe(i,1);
        end
    end

    % Nodal distributed forces
    P = zeros(nnodes,6);
    for i=1:length(PeS(:,1))
        P(PeS(i,2),PeS(i,3)) = PeS(i,1);
    end
    
    % Nodal body forces
    B = zeros(nnodes,6);
    for i=1:length(Be(:,1))
        B(Be(i,2),Be(i,3)) = Be(i,1);
    end

    % Assembly of all forces types
    b = zeros(4*6,1);
    p = zeros(4*6,1);
    fe = zeros(4*6,1);
    for e=1:nel
        b(:,e) = [B(Tn(e,1),:) B(Tn(e,2),:) B(Tn(e,3),:) B(Tn(e,4),:)]';
        p(:,e) = [P(Tn(e,1),:) P(Tn(e,2),:) P(Tn(e,3),:) P(Tn(e,4),:)]';
        fe(:,e) = Me(:,:,e)*b(:,e);
        for j=1:4
            % Element force vector computation
            fe(:,e) = fe(:,e) + S_4(e,j)*R(:,:,e)'*N_mass(:,:,e,j)'*N_mass(:,:,e,j)*R(:,:,e)*p(:,e);
        end
        
        % Global force vector assembly
        Idof = zeros(24,1);
        for j=1:6
            Idof(j,1) = 6*(Tn(e,1)-1)+j;
            Idof(6+j,1) = 6*(Tn(e,2)-1)+j;
            Idof(12+j,1) = 6*(Tn(e,3)-1)+j;
            Idof(18+j,1) = 6*(Tn(e,4)-1)+j;       
        end
        F(Idof,1) = F(Idof,1) + fe(:,e);
    end
end