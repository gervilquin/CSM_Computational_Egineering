function [u,FR] = SystemSolver(K,F,u,If,Ip)

    u(If,1) = K(If,If)\(F(If,1)-K(If,Ip)*u(Ip,1)); % Compute unknown displacements
    FR = K*u + F; % Compute reactions
    
end