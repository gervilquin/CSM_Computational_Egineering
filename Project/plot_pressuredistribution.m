function plot_pressuredistribution(X,n_l,n_u,p_inf,AoA)

    xu = unique(X(n_u(:,4),1));
    zu = unique(X(n_u(:,4),3));

    xl = unique(X(n_l(:,4),1));
    zl = unique(X(n_l(:,4),3));


    b = max(X(:,2)) - min(X(:,2));
    chord = max(X(:,1)) - min(X(:,1));

    y = 0*b;
    
    p_u = [];
    p_l = [];
    for i=1:length(xu)
        p_u(end+1) = ComputePressureSkin(AoA,xu(i),y,zu(i),chord,b,p_inf);
    end

    for i=1:length(xl)
        p_l(end+1) = ComputePressureSkin(AoA,xl(i),y,zl(i),chord,b,p_inf);
    end

    figure()
    hold on
    title("Pressure distribution at y = 0 ",'Interpreter','latex')
    plot(xu/chord,p_u)
    plot(xl/chord,p_l)
    plot(xl/chord,zeros(1,length(xl)),'color','k','HandleVisibility','off')
    legend(["Suction side","Pressure side"])
    grid on
    grid minor
    ylabel("Pressure [Pa]",'Interpreter','latex')
    xlabel("x/c","Interpreter","latex")
    hold off




end