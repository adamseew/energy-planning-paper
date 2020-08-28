
function [] = build_gdn(E, ke, varphi, x, y)
   
    dvarphi = gradient(varphi, [x y]);
    
    dpd = E * dvarphi - ke * varphi * dvarphi;
    
    [xp yp] = meshgrid(0:.5:25, 0:.5:25);

    points = reshape(cat(2, xp', yp'), [], 2);

    clear xp yp;

    figure;

    
    title('vector field with constant v');

    varphi_plot = ezplot(varphi, [0 25]);
    set(varphi_plot, 'Color', 'black', 'LineWidth', 2)

    xlim([0 25]);
    ylim([0 25]);

    hold on;

    dpdx = [];
    dpdy = [];
    
    for i = 1:size(points, 1)
    
        p = transpose(points(i, :));
        
        x = p(1);
        y = p(2);
        
        pd = double(subs(dpd));
    
        angle = atand(pd(2) / pd(1));
    
        if (isnan(angle)) % center of the ellipse
            continue;
        end
    
        %ha = annotation('arrow');
        %ha.Parent = gca;
        %ha.Color = [170 170 170] / 255;
    
        if (pd(1) < 0 && pd(2) < 0)
            angle = 180 + angle;
        end
    
        if (pd(1) < 0 && pd(2) > 0)
            angle = 180 + angle;
        end
  
        %ha.X = [p(1) p(1)+0.45*cosd(angle)];
        %ha.Y = [p(2) p(2)+0.45*sind(angle)];
        
        dpdx = [dpdx [p(1); p(1) + 0.45 * cosd(angle)]];
        dpdy = [dpdy [p(2); p(2) + 0.45 * sind(angle)]];
        

        %plot([p(1); p(1)+0.45*cosd(angle)],[p(2); p(2)+0.45*sind(angle)]); 
        
        %ha.LineWidth  = 1;
        %ha.HeadWidth  = 2;
        %ha.HeadLength = 4;
    end
    
    plot(dpdx, dpdy); 
    
end
