function [sigVM] = ComputeVonMissesStresses(Tn, u, Bs, Bmt, Bmn, Bb, R, nu, E, h)
    % Initialisation
    Idof = zeros(4*6,1);
    eps_b = zeros(3,size(Tn,1),4);
    eps_m = zeros(3,size(Tn,1),4);
    eps_s = zeros(2,size(Tn,1),4);
    sig_b = zeros(3,size(Tn,1),4);
    sig_m = zeros(3,size(Tn,1),4);
    sig_s = zeros(2,size(Tn,1),4);
    sigVM = zeros(size(Tn,1),4);

    for e=1:size(Tn,1)
        % Get strain components
        for j=1:6 % #ndofs
            Idof(j,1) = 6*(Tn(e,1)-1) + j;
            Idof(6+j,1) = 6*(Tn(e,2)-1) + j;
            Idof(12+j,1) = 6*(Tn(e,3)-1) + j;
            Idof(18+j,1) = 6*(Tn(e,4)-1) + j;
        end
        
        for k=1:4 % #ngauss
            eps_b(:,e,k) = Bb(:,:,e,k)*R(:,:,e)*u(Idof,1);
            eps_m(1:2,e,k) = Bmn(:,:,e)*R(:,:,e)*u(Idof,1);
            eps_m(3,e,k) = Bmt(:,:,e)*R(:,:,e)*u(Idof,1);
            eps_s(:,e,k) = Bs(:,:,e)*R(:,:,e)*u(Idof,1);
        end
        % Get stress components
        Cp = (E/(1-nu))*[ 1   nu      0   ;
                          nu  1       0   ;
                          0   0   (1-nu)/2];
        Cs = (E/(2*(1+nu)))*eye(2);

        for k=1:4
            sig_m(:,e,k) = Cp*eps_m(:,e,k);
            sig_s(:,e,k) = Cs*eps_s(:,e,k);
            sig_b(:,e,k) = Cp*h*0.5*eps_b(:,e,k);
            sig_plus = [sig_m(:,e,k) + sig_b(:,e,k); sig_s(:,e,k)]';
            sig_VMplus = (sig_plus(1)^2 + sig_plus(2)^2 - sig_plus(1)*sig_plus(2) + 3*(sig_plus(3)+sig_plus(4)+sig_plus(5)))^(0.5);
            sig_minus = [sig_m(:,e,k) - sig_b(:,e,k); sig_s(:,e,k)]';
            sig_VMminus = (sig_minus(1)^2 + sig_minus(2)^2 - sig_minus(1)*sig_minus(2) + 3*(sig_minus(3)+sig_minus(4)+sig_minus(5)))^(0.5);
            sigVM(e,k) = max(sig_VMplus,sig_VMminus);
        end
    end
end