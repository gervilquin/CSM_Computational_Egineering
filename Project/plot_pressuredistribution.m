function plot_pressuredistribution(X,n_l,n_u,p_inf,AoA)

    xu = unique(X(n_u(:,4),1));
    zu = unique(X(n_u(:,4),3));

    xl = unique(X(n_l(:,4),1));
    zl = unique(X(n_l(:,4),3));


    b = max(X(:,2)) - min(X(:,2));
    chord = max(X(:,1)) - min(X(:,1));

    y = 0*b;
    
    p_u = zeros(length(xl),2);
    p_l = zeros(length(xl),2);
    for i=1:length(xu)
        p_u(i,1) = ComputePressureSkin(AoA,xu(i),y,zu(i),chord,b,p_inf);
    end

    for i=1:length(xl)
        p_l(i,1) = ComputePressureSkin(AoA,xl(i),y,zl(i),chord,b,p_inf);
    end

    y = 0.9*b;
    for i=1:length(xu)
        p_u(i,2) = ComputePressureSkin(AoA,xu(i),y,zu(i),chord,b,p_inf);
    end

    for i=1:length(xl)
        p_l(i,2) = ComputePressureSkin(AoA,xl(i),y,zl(i),chord,b,p_inf);
    end


    figure()
    hold on
    title("Pressure distribution at y = 0 m and y = 2.7 m  ",'Interpreter','latex')
    plot(xu/chord,p_u(:,1),'DisplayName','Suction side y = 0')
    plot(xl/chord,p_l(:,1),'DisplayName','Pressure side y = 0')
    plot(xu/chord,p_u(:,2),'DisplayName','Suction side y = 2.7m')
    plot(xl/chord,p_l(:,2),'DisplayName','Pressure side y = 2.7m')
    plot(xl/chord,zeros(1,length(xl)),'color','k','HandleVisibility','off')
    %legend(["Suction side","Pressure side"])
    legend('Interpreter','latex')
    set(gca, 'YDir','reverse')
    grid on
    grid minor
    ylabel("Pressure [Pa]",'Interpreter','latex')
    xlabel("x/c","Interpreter","latex")
    hold off




end