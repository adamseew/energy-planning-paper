
E = [ 0 1; -1 0 ];

ke = 0.1;

varphi = @(p) ([1 0]*p - 5)^2/3^2 + ([0 1]*p - 5)^2/4^2 - 1;

dvarphi = @(p) [2*[1 0]*p/3^2 - 10/3^2; 2*[0 1]*p/4^2 - 10/4^2];

dpd = @(p, ke) E * dvarphi(p) - ke * varphi(p) * dvarphi(p);

[x y] = meshgrid(0:.5:10, 0:.5:10);

points = reshape(cat(2, x', y'), [], 2);

clear x y;

figure;

varphi_plot = ezplot('(x-5)^2/3^2+(y-5)^2/4^2=1',[0 10]);
set(varphi_plot, 'Color', 'black', 'LineWidth', 2)


xlim([0 10]);
ylim([0 10]);

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
    ha.HeadWidth  = 3;
    ha.HeadLength = 5;
end




