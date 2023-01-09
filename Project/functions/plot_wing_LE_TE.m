function plot_wing_LE_TE(X,u,u2,u6,I_LE,I_TE)

    DOF = 6;
    u_matrix = vector2matrix(u,DOF);
    u2_matrix = vector2matrix(u2,DOF);
    u6_matrix = vector2matrix(u6,DOF);




    figure()
    title("Leading Edge",'Interpreter','latex')
    hold on
    plot(X(I_LE,2),u_matrix(I_LE,3),'-x','color','black','DisplayName','inverse')
    plot(X(I_LE,2),u2_matrix(I_LE,3),'color','blue','DisplayName','2 modes')
    plot(X(I_LE,2),u6_matrix(I_LE,3),'color','green','DisplayName','6 modes')
    legend("Interpreter","latex",'Location','northwest')
    grid on
    grid minor
    xlabel("Y[m]",'Interpreter','latex')
    ylabel("$u_y[m]$",'Interpreter','latex')
    hold off

    figure()
    title("Trailing Edge",'Interpreter','latex')
    hold on
    plot(X(I_TE,2),u_matrix(I_TE,3),'-x','color','black','DisplayName','inverse')
    plot(X(I_TE,2),u2_matrix(I_TE,3),'color','blue','DisplayName','2 modes')
    plot(X(I_TE,2),u6_matrix(I_TE,3),'color','green','DisplayName','6 modes')
    legend("Interpreter","latex",'Location','northwest')
    grid on
    grid minor
    xlabel("Y[m]",'Interpreter','latex')
    ylabel("$u_y[m]$",'Interpreter','latex')
    hold off

    figure()
    title("Leading and Trailing Edges",'Interpreter','latex')
    hold on
    plot(X(I_LE,2),u_matrix(I_LE,3),'-x','color','black','DisplayName','LE inverse')
    plot(X(I_LE,2),u2_matrix(I_LE,3),'color','blue','DisplayName','LE 2 modes')
    plot(X(I_LE,2),u6_matrix(I_LE,3),'color','green','DisplayName','LE 6 modes')

    plot(X(I_TE,2),u_matrix(I_TE,3),'--o','color','red','DisplayName','TE inverse')
    plot(X(I_TE,2),u2_matrix(I_TE,3),'--','color','magenta','DisplayName','TE 2 modes')
    plot(X(I_TE,2),u6_matrix(I_TE,3),'--','color','cyan','DisplayName','TE 6 modes')
    legend("Interpreter","latex",'Location','northwest')
    grid on
    grid minor
    xlabel("Y[m]",'Interpreter','latex')
    ylabel("$u_y[m]$",'Interpreter','latex')
    hold off


end