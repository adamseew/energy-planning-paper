%OP1
% Simulation of the Energy-Aware Dynamic Mission Planning Algorithm
%
% Vector field simulation from a defined TEE, plot and possible animation


%% Build the model 

fprintf(['[ OP1 ] Vector field simulation from a defined TEE, plot and\n' ...
         '        possible animation\n']);
disp( '[     ] Build the model');

syms x y;

varphi = input('[   ? ] Type the TEE: ');

if isempty(varphi)
   
    min = 0;
    max = 400;
    varphi = ((x - 200) * cosd(12) + (y - 200) * sind(12))^2 / 120^2 + ...
        ((x - 200)*sind(12) - (y - 200) * cosd(12))^2 / 170^2 - 1; 
    % default TEE, an ellipse
else
    
    min = input('[   ? ] The min value: ');
    max = input('[   ? ] And the max value: ');
end

% create all the points used to build the vector field

[xp yp] = meshgrid(...
    min:(abs(min) + abs(max)) / 50:max, min:(abs(min) + abs(max)) / 50:max...
                  );
              
points = reshape(cat(2, xp', yp'), [], 2);
clear xp yp;

E = [0 1; -1 0]; % default rotation matrix

ke = input('[   ? ] Gain to adjust speed of convergence: ');
if isempty(ke)
   ke = 1.5; % default gain to adjust speed of convergence
end

% Builds the vector field
[dpd, pdangle] = build_gdn2(E, ke, varphi, points);


%% Plot and animate

disp('[     ] Plot and animate');
fprintf(['[   ! ] Dependencies: ration matrix E, gain ke, TEE varphi,\n' ...
        '        points, the angle of the vector field pdangle, the min\n' ...
        '        and max value\n']);

if strcmp(questdlg('Show the animation?', 'Animate', 'Yes', 'No', 'Yes'),...
           'Yes'...
         )

    answer = inputdlg({'start coordinate x:', 'start coordinate y:'}, ...
                      'Starting point', [1 35], {'0', '0'} ...
                     );

    if isempty(answer)
        strp = [0; 0];
    else
        strp = str2double(answer);
    end
    
    clear answer;
    
    andu = input('[   ? ] Animation duration: ');
    if isempty(andu)
        andu = 10000; % default animation duration
    end

    s = input('[   ? ] Cruise speed [m/s]: ');
    if isempty(s)
        s = 20; % default cruise speed
    end

    plot_gdn2(E, ke, varphi, points, pdangle, ...
        min, max, andu, s, strp);
else
    plot_gdn2(E, ke, varphi, points, pdangle, ...
        min, max, 0, 0, [0; 0]);
end

