function Up = AssignPrescribedDisplacements(y_fix_nodes,X)
    DOF = 6;
    index_nodes_yp = find(X(:,2)==y_fix_nodes);

    Up =zeros(length(index_nodes_yp)*DOF,3);

    for i =1:length(index_nodes_yp)
        index_x_up = (i-1)*DOF;

                        %       u   n                   j 
        Up(index_x_up+1,:) = [  0   index_nodes_yp(i)   1];
        Up(index_x_up+2,:) = [  0   index_nodes_yp(i)   2];
        Up(index_x_up+3,:) = [  0   index_nodes_yp(i)   3];
        Up(index_x_up+4,:) = [  0   index_nodes_yp(i)   4];
        Up(index_x_up+5,:) = [  0   index_nodes_yp(i)   5];
        Up(index_x_up+6,:) = [  0   index_nodes_yp(i)   6];
    end

end