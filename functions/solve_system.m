function [u,FR] = solve_system(u,K,F,If,Ip)
    
    % displacements and rotations at free DOFs
    u(If,1) = K(If,If)\(F(If,1) - K(If,Ip)*u(Ip,1));

    % reaction forces adn moments at DOFs
    FR = K*u + F;

end