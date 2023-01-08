function [G,j_p,A,J,Iy,Iz,ky,kz,kt] = beam_section_properties(E,nu,d)
    
    %Shear modulus of isotropic material
    G = E/(2*(1+nu));

    % Orientation of the y axis
    j_p = [0,0,1];

    % Area section
    A = pi/4 * d^2;

    % Inertia
    % ^ z
    % |
    % |
    % |--------> y


    Iy = A*(d/2)^2/4;
    Iz = Iy;
    J = Iy + Iz;

    % Shear torsion correction factors -> Imposed
    ky = 6/7;
    kz = ky;
    kt = 1;






end