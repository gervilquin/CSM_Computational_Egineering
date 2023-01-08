function [u,R] = SystemSolver(K,F,u,If,Ip)

    u(If,1) = K(If,If)\(F(If,1)-K(If,Ip)*u(Ip,1)); % Compute unknown displacements
    R = K*u + F; % Compute reactions
    
end