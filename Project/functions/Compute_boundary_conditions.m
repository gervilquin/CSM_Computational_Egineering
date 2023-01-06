function [u,If,Ip] = Compute_boundary_conditions(Ne,Up)
    
    % Initialization
    Ndof = (Ne+1)*6;
    u = zeros(Ndof,1);
    Ip = zeros(length(Up(:,1)),1);

    % Prescribed and free DOFs
    for p = 1:length(Up(:,1))
        %disp(" >")
        Ip(p) = 6*(Up(p,2)-1) + Up(p,3);
        u(Ip(p),1) = Up(p,1);
    end
    If = setdiff(1:Ndof,Ip);


end