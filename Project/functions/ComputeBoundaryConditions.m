function [If, Ip, u] = ComputeBoundaryConditions(nDOFs,Up)
    % Initialization
    u = zeros(nDOFs,1);
    Ip = zeros(length(Up),1); % DOFs with prescribed displacement

    for i=1:length(Up)
        Ip(i) = 6*(Up(i,2)-1) + Up(i,3);
        u(Ip(i),1) = Up(i,1);
    end
    If = setdiff(1:nDOFs,Ip); % DOFs with free displacement
end