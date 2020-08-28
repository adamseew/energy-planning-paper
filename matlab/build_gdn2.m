
function [] = build_gdn(E, ke, varphi, min, max)
    
    syms x y;
    
    dvarphi = gradient(varphi, [x y]);
    
    dpd = E * dvarphi - ke * varphi * dvarphi;
    
    [xp yp] = meshgrid(min:(abs(min) + abs(max)) / 50:max, min:(abs(min) + abs(max)) / 50:max);

    points = reshape(cat(2, xp', yp'), [], 2);

    clear xp yp;

    figure;
    
    title('vector field with constant v');

    varphi_plot2 = fimplicit(varphi);
    set(varphi_plot2, 'Color', [.85 .85 .85], 'LineWidth', 10);
    
    hold on;
    
    varphi_plot = fimplicit(varphi);
    set(varphi_plot, 'Color', 'black', 'LineWidth', 1.2);
    
    xlim([min max]);
    ylim([min max]);
        
    x = transpose(points(:, 1));
    y = transpose(points(:, 2));
        
    pd = transpose(double(subs(dpd)));
    
    angle = atand(pd(:,2) ./ pd(:,1));
    
    angle(and(pd(:,1) < 0, pd(:,2) < 0)) = angle(and(pd(:,1) < 0, pd(:,2) < 0)) + 180;
   
    %if (pd(1) < 0 && pd(2) < 0)
    
    angle(and(pd(:,1) < 0, pd(:,2) > 0)) = angle(and(pd(:,1) < 0, pd(:,2) > 0)) + 180;
    
    %if (pd(1) < 0 && pd(2) > 0)
     
    angle = transpose(angle);
    
    % just plots
    
    coeff = 0.45 * (max - min) / 25;
    coeffm = 0.20 * (max - min) / 25;
    coeffM = 0.30 * (max - min) / 25;
    
    dpdx = [x; x + coeff * cosd(angle)];
    dpdy = [y; y + coeff * sind(angle)];
    
    arrowlxx = [x + coeffM * cosd(angle); dpdx(2, :)];
    arrowlxy = [y + coeffm * sind(angle); dpdy(2, :)];
    
    arrowrxx = [x + coeffm * cosd(angle); dpdx(2, :)];
    arrowrxy = [y + coeffM * sind(angle); dpdy(2, :)];
    
    plot(dpdx, dpdy, 'Color', [0.3010, 0.7450, 0.9330]);
    
    plot(arrowlxx, arrowlxy, 'Color', [0.3010, 0.7450, 0.9330]);
    plot(arrowrxx, arrowrxy, 'Color', [0.3010, 0.7450, 0.9330]);
    
    % animation of 100 s
    
    set(gcf, 'Position', get(0, 'Screensize'));
    n = 1000000;
    xd = 400;
    yd = 400;
    s = 20;
    xdd = xd;
    ydd = yd;
    for i=1:n
        x = xd;
        y = yd;
        pd = transpose(double(subs(dpd)));
        
        angle = atand(pd(:,2) ./ pd(:,1));
        angle(and(pd(:,1) < 0, pd(:,2) < 0)) = angle(and(pd(:,1) < 0, pd(:,2) < 0)) + 180;
        angle(and(pd(:,1) < 0, pd(:,2) > 0)) = angle(and(pd(:,1) < 0, pd(:,2) > 0)) + 180;
    
        angle = transpose(angle);
        
        xd = [xd + s * .1 * cosd(angle)];
        yd = [yd + s * .1 * sind(angle)];
        
        xdd = [xdd; xd];
        ydd = [ydd; yd];
        
        plot(xdd,ydd,'Color','r','LineWidth', 1.2)
        
        pause(.1)
    end
    
end
