
E = [ 0 1; -1 0 ];

ke = 1.5;

varphi = @(p) ([1 0]*p - 12.5)^2/6^2 + ([0 1]*p - 12.5)^2/8^2 - 1;

dvarphi = @(p) [2*[1 0]*p/6^2 - 25/6^2; 2*[0 1]*p/8^2 - 25/8^2];

dpd = @(p, ke) E * dvarphi(p) - ke * varphi(p) * dvarphi(p);

[x y] = meshgrid(0:.5:25, 0:.5:25);

points = reshape(cat(2, x', y'), [], 2);

clear x y;

figure;

title('vector field (aka desired velocity vector)');

xlim([0 25]);
ylim([0 25]);

hold on;

for i = 1:size(points, 1)
    
    p = transpose(points(i, :));
    
    pd = dpd(p, ke);
    
    ha = annotation('arrow');
    ha.Parent = gca;
    ha.Color = [170 170 170] / 255;
    
    ha.X = [p(1) p(1)+pd(1)];
    ha.Y = [p(2) p(2)+pd(2)]; 

    ha.LineWidth  = 1;
    ha.HeadWidth  = 2;
    ha.HeadLength = 4;
end

figure;

title('vector field with constant v');

varphi_plot = ezplot('(x-12.5)^2/6^2+(y-12.5)^2/8^2=1',[0 25]);
set(varphi_plot, 'Color', 'black', 'LineWidth', 2)

[x y] = meshgrid(0:0.5:25, 0:0.5:25);

points = reshape(cat(2, x', y'), [], 2);

xlim([0 25]);
ylim([0 25]);

hold on;

for i = 1:size(points, 1)
    
    p = transpose(points(i, :));
    
    pd = dpd(p, ke);
    
    angle = atand(pd(2) / pd(1));
    
    if (isnan(angle)) % center of the ellipse
        continue;
    end
    
    ha = annotation('arrow');
    ha.Parent = gca;
    ha.Color = [170 170 170] / 255;
    
    if (pd(1) < 0 && pd(2) < 0)
        angle = 180 + angle;
    end
    
    if (pd(1) < 0 && pd(2) > 0)
        angle = 180 + angle;
    end
  
    ha.X = [p(1) p(1)+0.45*cosd(angle)];
    ha.Y = [p(2) p(2)+0.45*sind(angle)];

    ha.LineWidth  = 1;
    ha.HeadWidth  = 2;
    ha.HeadLength = 4;
end





