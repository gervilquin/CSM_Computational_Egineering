function [F] = ComputeFvectorBeam(nDOFs,Tn,Fe,PeB,Be,Me,R,w,Nek,le,FeNotComputed)
   %Initialization
    nel = size(Tn,1);
    nnodes = 4*nel;
    F = zeros(nDOFs,1);

    % Point loads
    if FeNotComputed == true 
        for i= 1:length(Fe(:,1))
            F(6*(Fe(i,2)-1) + Fe(i,3),1) = F(6*(Fe(i,2)-1) + Fe(i,3),1) + Fe(i,1);
        end
    end

    % Nodal distributed loads
    P = zeros(nnodes,6);
    for i=1:length(PeB(:,1))
        P(PeB(i,2),PeB(i,3)) = PeB(i,1);
    end

    % Nodal body forces
    B = zeros(nnodes,6);
    for i=1:length(Be(:,1))
        B(Be(i,2),Be(i,3)) = Be(i,1);
    end

    % Assembly of all forces types
    b = zeros(2*6,1);
    p = zeros(2*6,1);
    fe = zeros(2*6,1);
    for e=1:nel
        b(:,e) = [B(Tn(e,1),:),B(Tn(e,2),:)]';
        p(:,e) = [P(Tn(e,1),:),P(Tn(e,2),:)]';
        fe(:,e) = Me(:,:,e)*b(:,e);
        for j=1:2
            % Element force vector computation
            fe(:,e) = fe(:,e) + w(j)*le(e)*R(:,:,e)'*Nek(:,:,e,j)'*Nek(:,:,e,j)*R(:,:,e)*p(:,e)/2;
        end

        % Global force vector assembly
        Idof = zeros(12,1);
        for j =1:6
            Idof(j) = 6*(Tn(e,1)-1) + j;
            Idof(6+j) = 6*(Tn(e,2)-1) + j;
        end
        F(Idof,1) = F(Idof,1) + fe(:,e); 
    end
end