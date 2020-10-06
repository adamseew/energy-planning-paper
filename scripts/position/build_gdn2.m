function [dpd, pdangle] = build_gdn2(E, ke, varphi, points)
%BUILD_GDN2 Creates the vector field
%
% Inputs:
%   E:      rotation matrix
%   ke:     gain to adjust speed of convergence
%   varphi: TEE
%   points: vector of points where the function generates the vector field
%   
%
% Outputs:
%   dpd:    desired velocity vector (generated from points)
%   pdangle:angle saying where the desired velocity vector is pointing
%
    syms x y;
    
    dvarphi = gradient(varphi, [x y]); % gradient
    
    dpd = E * dvarphi - ke * varphi * dvarphi;
    
    x = transpose(points(:, 1));
    y = transpose(points(:, 2));
    
    pd = transpose(double(subs(dpd)));
    
    angle = atand(pd(:, 2) ./ pd(:, 1));
    
    % One has to do the following operation, as the sign is not preserved
    % in the atand function above
    
    % Corresponds to if (pd(1) < 0 && pd(2) < 0) ...
    angle(and(pd(:, 1) < 0, pd(:, 2) < 0)) = ...
        angle(and(pd(:, 1) < 0, pd(:, 2) < 0)) + 180;
    
    % And this corresponds to if (pd(1) < 0 && pd(2) > 0) ...
    angle(and(pd(:, 1) < 0, pd(:, 2) > 0)) = ...
        angle(and(pd(:, 1) < 0, pd(:, 2) > 0)) + 180;
     
    angle = transpose(angle);
    pdangle = angle;
    dpd = pd;
    
end

