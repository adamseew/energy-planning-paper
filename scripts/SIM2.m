
% simulation with fixed parameters (controls)
% the simulations are static and they generate the plots of
%   - path (x,y)
%   - energy signal (time, power) [h in Equation (2)]
%   - estimated energy output (same as before) [y in Equation (3)]

% asking data about the path

answer = inputdlg(...
    {'vehicle speed [m/s]:', 'wind speed [m/s]:', ...
     'vehicle direction [angles]:', 'wind direction [angles]:', ...
     'start coordinate x [m]:', 'start coordinate y [m]:', ...
     'max power [W]:', 'min power [W]:', 'triggering point radius [m]' ...
    }, ...
    'path initialization', [1 40], ...
    {'20', '5', '270', '90', '-100', '220', '60', '30', '10'}); % asking for initial data

if isempty(answer)
    strp = [20; 5; 270; 90; -100; 220; 60; 30; 10]; % default initial data
else
    strp = str2double(answer);
end

clear answer;

% storing data

vehspeed = strp(1); 
windspeed = strp(2); 
vehdir = strp(3); 
winddir = strp(4);
startx = strp(5);
starty = strp(6);
maxpw = strp(7);
minpw = strp(8); 
trigeps = strp(9);

clear strp;

% asking data about the algorithm

answer = inputdlg(...
    {'order r:', 'epsilon:', 'horizon N:'}, ...
    'model and algorithm initialization', [1 40], ...
    {'3', '1', '10'}); % asking for initial data

if isempty(answer)
    strp = [3; 10; 1; 10]; % default initial data
else
    strp = str2double(answer);
end

clear answer;

% stroing data

r = strp(1);
xi = strp(2);
eps = strp(3);
N = strp(4);

clear strp;

% asking for the plan specifization 

disp('plan specification');
[bn, folder] = uigetfile('.pln');
if bn == 0
    fid = fopen('../data/simulation3/plan_specification.pln', 'rt');
else
    fid = fopen(fullfile([folder, bn]), 'rt');
end

path = []; % trajectories

triggers = []; % triggering points

i = -1;
while true
    thisline = fgetl(fid);
    if ~ischar(thisline) 
        break; 
    elseif isempty(thisline)
        continue;
    end  % end of file
    
    if contains(thisline, "[")
        thisline(1:1) = [];
        thisline(end:end) = [];
        thisline = split(thisline, ',');
        % in this simulation parameters are ignored
    elseif contains(thisline, ",")
        thisline = split(thisline, ",");
        triggers = [triggers; str2double(thisline(1)) str2double(thisline(2))];
    else
        path = [path; convertCharsToStrings(thisline)];
        i = i + 1;
    end
end

fclose(fid);

% building the model

[A C] = build_model(1, r);
syms q [2*r+1 1];
delta = .01;

% discretizing 

Ad = A*delta+eye(2*r+1);

% iterating trajectories in the plan to get the period

for traj = transpose(path) % per each trajectory
        
    traj = split(traj, ";");
    traj = str2sym(traj(3));
    
    
end

% model building function

function [A, C] = build_model(omega, r)

    m = 2*r+1;

    Aj = @(omega, j) [0 omega*j ; -omega*j 0];
    A  = zeros(m);
    A(1,1) = 0;
    C = [1];

    for i = 1:r
        A(2*i : 2*i + 1, 2*i : 2*i + 1) = Aj(2*pi/1, i); 
        C = [C 1 0];
    end
    
    C = 1/1 * C;

end
