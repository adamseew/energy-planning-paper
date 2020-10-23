function [dpdx, dpdy] = plot_gdn2(E, ke, varphi, points, pdangle, min, max, andu, s, p)
% PLOT_GDN2 Plots the vector field, and eventually animates it
%
% Inputs:
%   E:      rotation matrix
%   ke:     gain to adjust speed of convergence
%   varphi: TEE
%   points: vector of points where the function generates the vector field
%   pdangle:vector of desired velocity angles
%   min:    position (x and y) of the start of the vector field plot
%   max:    position of the end of the vector field plot
%   andu:   animation duration (use 0 if no animation required)
%   s:      cruise speed in m/s
%   p:      starting point (for the anumation)
%
    
    % Just a nice effect nothing more...
    varphi_plot2 = fimplicit(varphi);
    set(varphi_plot2, 'Color', [.85 .85 .85], 'LineWidth', 10);
    
    hold on;
    
    varphi_plot = fimplicit(varphi);
    set(varphi_plot, 'Color', 'black', 'LineWidth', 1.2);
    
    xlim([min max]);
    ylim([min max]);
        
    x = transpose(points(:, 1));
    y = transpose(points(:, 2));
        
    
    % And just a crazy workaround to plot the arrows...
    coeff = 0.45 * (max - min) / 25;
    coeffm = 0.20 * (max - min) / 25;
    coeffM = 0.30 * (max - min) / 25;
    
    dpdx = [x; x + coeff * cosd(pdangle)];
    dpdy = [y; y + coeff * sind(pdangle)];
    
    arrowlxx = [x + coeffM * cosd(pdangle); dpdx(2, :)];
    arrowlxy = [y + coeffm * sind(pdangle); dpdy(2, :)];
    
    arrowrxx = [x + coeffm * cosd(pdangle); dpdx(2, :)];
    arrowrxy = [y + coeffM * sind(pdangle); dpdy(2, :)];
    
    plot(dpdx, dpdy, 'Color', [0.3010, 0.7450, 0.9330]);
    
    plot(arrowlxx, arrowlxy, 'Color', [0.3010, 0.7450, 0.9330]);
    plot(arrowrxx, arrowrxy, 'Color', [0.3010, 0.7450, 0.9330]);
    
    % Here is the eventual animation
    xd = [1 0] * p;
    yd = [0 1] * p;
    
    xdd = xd;
    ydd = yd;
    
    for i=1:andu
        x = xd;
        y = yd;
        [dpd, pdangle] = build_gdn2(E, ke, varphi, [x y]);
        pd = transpose(double(subs(dpd)));
        
        pdangle = transpose(pdangle);
        
        xd = [xd + s * .1 * cosd(pdangle)];
        yd = [yd + s * .1 * sind(pdangle)];
        
        xdd = [xdd; xd];
        ydd = [ydd; yd];
        
        plot(xdd, ydd, 'Color', 'r', 'LineWidth', 1.2)
        
        pause(.1)
    end
end

