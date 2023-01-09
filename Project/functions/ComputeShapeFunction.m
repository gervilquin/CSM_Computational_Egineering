function [N,derN,wk] = ComputeShapeFunction(nG)
% Nodal coefficients determine which node shape function is being used. The
% position of the Gauss points determine which points of those functions
% are we using to asses the integral.
    N = zeros(4,nG); 
    derN = zeros(4,2,nG); % j_dim -> 2 derivation variables (1:xi, 2:eta)
    a_nodal_coef = [-1 1 1 -1]; % Coefficients for nodes 1,2,3,4.
    b_nodal_coef = [-1 -1 1 1]; % idem

    if nG == 1
        wk = 4;
        GP_local_coord = [0 0];
    elseif nG == 4
        wk = [1 1 1 1]';
        GP_local_coord = [ -1/sqrt(3) -1/sqrt(3);
                           +1/sqrt(3) -1/sqrt(3);
                           +1/sqrt(3) +1/sqrt(3);
                           -1/sqrt(3) +1/sqrt(3)]; % [xi, eta]
    end

    for i=1:4
        for j=1:nG
            N(i,j) = (1/4)*(1+a_nodal_coef(i)*GP_local_coord(j,1))*(1+b_nodal_coef(i)*GP_local_coord(j,2));
            derN(i,1,j) = (1/4)*a_nodal_coef(i)*(1+b_nodal_coef(i)*GP_local_coord(j,2));
            derN(i,2,j) = (1/4)*b_nodal_coef(i)*(1+a_nodal_coef(i)*GP_local_coord(j,1));
        end
    end
end