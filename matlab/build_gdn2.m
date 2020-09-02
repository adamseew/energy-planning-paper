
function [dpd, pdangle] = build_gdn2(E, ke, varphi, points)
    
    syms x y;
    
    % gradient
    dvarphi = gradient(varphi, [x y]);
    
    dpd = E * dvarphi - ke * varphi * dvarphi;
    
    x = transpose(points(:, 1));
    y = transpose(points(:, 2));
    
    pd = transpose(double(subs(dpd)));
    
    angle = atand(pd(:, 2) ./ pd(:, 1));
    
    angle(and(pd(:, 1) < 0, pd(:, 2) < 0)) = angle(and(pd(:, 1) < 0, pd(:, 2) < 0)) + 180;
    
    % corresponds to if (pd(1) < 0 && pd(2) < 0) ...
    
    angle(and(pd(:, 1) < 0, pd(:, 2) > 0)) = angle(and(pd(:, 1) < 0, pd(:, 2) > 0)) + 180;
    
    % corresponds to if (pd(1) < 0 && pd(2) > 0) ...
     
    angle = transpose(angle);
    pdangle = angle;
    dpd = pd;
end
