function [Sa,Ss,St,Sb,Fx,Fy,Fz,Mx,My,Mz] = Compute_interal_forces_strain(Ne,Tn,u,Ba,Bs,Bt,Bb,R,Ka,Kb,Ks,Kt)

    % Initialization
    Sa = zeros(1,Ne);
    Ss = zeros(2,Ne);
    St = zeros(1,Ne);
    Sb = zeros(2,Ne);

    Fx = zeros(1,Ne); 
    Fy = zeros(1,Ne); 
    Fz = zeros(1,Ne); 

    Mx = zeros(1,Ne); 
    My = zeros(1,Ne); 
    Mz = zeros(1,Ne); 

    for e = 1:length(Tn)
        % Get element displacements
        Idof = zeros(12,1);
        for j =1:6
            Idof(j,1) = 6*(Tn(e,1)-1)+j;
            Idof(6+j,1) = 6*(Tn(e,2)-1)+j;
        end
        ue = u(Idof);

        % Get each strain components
        Sa(1,e) = Ba(1,:,e)*R(:,:,e)*ue;
        Ss(:,e) = Bs(:,:,e)*R(:,:,e)*ue;
        St(1,e) = Bt(1,:,e)*R(:,:,e)*ue;
        Sb(:,e) = Bb(:,:,e)*R(:,:,e)*ue;

        % Get internal forces and moments at each element node
        fint(:,e) = R(:,:,e)*(Ka(:,:,e)+Kb(:,:,e)+Ks(:,:,e)+Kt(:,:,e))*ue;
        
        Fx(:,e) = fint(1,e)-fint(7,e);
        Fy(:,e) = fint(2,e)-fint(8,e);
        Fz(:,e) = fint(3,e)-fint(9,e);

        Mx(:,e) = fint(4,e)-fint(10,e);
        My(:,e) = fint(5,e)-fint(11,e);
        Mz(:,e) = fint(6,e)-fint(12,e);


    end
    disp('A')

end