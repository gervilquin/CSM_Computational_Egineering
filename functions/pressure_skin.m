function p = pressure_skin(AoA,x,y,z,c,b,p_inf)
    
    if y <= 0.8*b
        P_y = p_inf;
    elseif y > 0.8*b
        P_y = p_inf*cos(pi/2 * (y-0.8*b)/(0.2*b));
    end

    if z >= 0
        p = -4*AoA*P_y*((1-x/c)^4 + sqrt(1-x/c));
    else
        p = -4*AoA*P_y*((1-x/c)^4 -1/4*sqrt(1-x/c));
    end
end 