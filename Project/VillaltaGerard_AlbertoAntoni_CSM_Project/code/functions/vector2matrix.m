function u = vector2matrix(u_in,DOF)
    
    Ndof = length(u_in);
    u = zeros(Ndof/DOF,DOF);

    for i = 1:DOF
        u(:,i) = u_in(i:DOF:Ndof-(DOF-i));
    end

end