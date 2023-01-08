function Pe = surface_forces(X,p_inf,n_u,n_l,AoA)
    
    N_skin_nodes = length(n_u(:,4))+length(n_l(:,4));
    chord = max(X(:,1));
    span = max(X(:,2));

    Pe = zeros(N_skin_nodes*3,3);

    

    for i = 1:length(n_l(:,1))
        p_node = pressure_skin(AoA,X(n_l(i,4),1),X(n_l(i,4),2),X(n_l(i,4),3),chord,span,p_inf);
        index_Pe = (i-1)*3;

                    %       u               n        j
        Pe(index_Pe+1,:) = [p_node*n_l(i,1) n_l(i,4) 1];
        Pe(index_Pe+2,:) = [p_node*n_l(i,2) n_l(i,4) 1];
        Pe(index_Pe+3,:) = [p_node*n_l(i,3) n_l(i,4) 1];
    end


    for i = 1:length(n_u(:,1))


        p_node = pressure_skin(AoA,X(n_u(i,4),1),X(n_u(i,4),1),X(n_u(i,4),1),chord,span,p_inf);
        index_Pe = (i-1)*3 + length(n_l(:,1))*3;

                    %       u               n        j
        Pe(index_Pe+1,:) = [p_node*n_u(i,1) n_u(i,4) 1];
        Pe(index_Pe+2,:) = [p_node*n_u(i,2) n_u(i,4) 1];
        Pe(index_Pe+3,:) = [p_node*n_u(i,3) n_u(i,4) 1];
    end




end