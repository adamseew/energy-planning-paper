
% simulation with fixed parameters (controls)
% the simulations are static and they generate the plots of
%   - path (x,y)
%   - energy signal (time, power) [h in Equation (2)]
%   - estimated energy output (same as before) [y in Equation (3)]


%% initialization

% asking data about the path

answer = inputdlg(...
    {'vehicle speed [m/s]:','wind speed [m/s]:',...
     'vehicle direction [angles]:','wind direction [angles]:',...
     'start coordinate x [m]:','start coordinate y [m]:',...
     'max power [W]:','min power [W]:','triggering point radius [m]'...
    }, ...
    'path initialization',[1 40],...
    {'20','5','270','90','-100','220','60','30','10'}); % asking for initial data

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
    {'order r:','epsilon:','horizon N:'},...
    'model and algorithm initialization',[1 40],...
    {'3','1','10'}); % asking for initial data

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
[bn folder] = uigetfile('.pln');
if bn == 0
    fid = fopen('../data/simulation3/plan_specification.pln','rt');
else
    fid = fopen(fullfile([folder bn]),'rt');
end


%% parsing plan

path = []; % trajectories

trigs = []; % triggering points

i = -1;
sp = [];
while true
    thisline = fgetl(fid);
    if ~ischar(thisline) 
        break; 
    elseif isempty(thisline)
        continue;
    end  % end of file
    
    if i == -1 % first line contains the shift
        thisline = split(thisline,',');
        sp = [str2double(thisline(1)); str2double(thisline(2))];
        i = i+1;
        continue;
    end
    
    if contains(thisline,"[")
        thisline(1:1) = [];
        thisline(end:end) = [];
        thisline = split(thisline,',');
        % in this simulation parameters are ignored
    elseif contains(thisline,",")
        thisline = split(thisline,",");
        trigs = [trigs; str2double(thisline(1)) str2double(thisline(2))];
    else
        path = [path; convertCharsToStrings(thisline)];
        i = i+1;
    end
end

fclose(fid);


%% building the model

[A C] = build_model(2*pi/1, r);
syms q [2*r+1 1];
delta = .01;

% discretizing 

Ad = A*delta+eye(2*r+1);


%% model guesses

q0 = ones(size(q,1),1);

P0 = eye(size(q,1));
    
% process noise and sensor noise

Q = P0;
R = 1;    


%% iterating trajectories in the plan to get the constant n 
% (to measure the period)

d = [];
p = [0; 0];
n = 0;

for traj = transpose(path) % per each trajectory
        
    traj = split(traj,";");
    traj = str2sym(traj(3));
    
    x = p(1);
    y = p(1);
    
    di = double(subs(traj));
    
    x = p(1) - sp(1);
    y = p(2) - sp(2);
    
    if ismember(double(subs(traj)),d)
        break;
    else
        d = [d di];
    end
    
    n = n + 1;
end


%% running the simulation
    
k = 0;
i = 1;
    
pos = [];
pow = [];
    
nowpos = [startx starty];
pdanglelist = [];
    
% contribution of the wind. This doesn't change as winspeed and
% direction is constant
windx = windspeed * cosd(winddir);
windy = windspeed * sind(winddir);
    
y = []; % model output
q = []; % model state

for traj = transpose(path)
        
    traj = split(traj, ";");
            
    if contains(traj(2), '1')
        E = [0 1; -1 0];
    else
        E = [0 -1; 1 0];
    end
        
    ke = str2double(traj(1));
        
    while true

        pos = [pos; nowpos];
            
        % vector field
        [dpd, pdangle] = build_gdn2(E, ke, str2sym(traj(3)), nowpos); 
                        
        pdanglelist = [pdanglelist; pdangle];
        
        % dpd is the ideal direction... We have wind though
        % now if the time is sampled every second, the below expression
        % actually indicates the offset from the original location
            
        posx = .1 * vehspeed * cosd(pdangle);
        posy = .1 * vehspeed * sind(pdangle);
        
        % new position
        nowpos = [nowpos(1) + windx + posx, nowpos(2) + windy + posy]; 
            
        % reached the triggering point, going to TEE i+1
        if all(abs(nowpos - trigs(i,:)) <= trigeps) 
            break;
        end
                 
        % produce a simulated energy value
        % first, let's calculate the angle between the two vectors
        angle = acosd(...
            (posx * windx + posy * windy) / ...
            (sqrt(posx^2 + posy^2) * sqrt(windx^2 + windy^2))...
                         );
        % angle is nan? no windspeed then; so the contribution of the
        % wind is null
        if isnan(angle)
            angle = 0;
        end
            
        % if the angle is 0, lowest energy (tail wind). 180, heighest
        % (head wind)
        simpow = minpw + (maxpw - minpw) * .5 * abs(1 - cosd(angle)) * ... 
            windspeed / vehspeed; % winddir and windspeed effect
            
        pow = [pow; simpow];
                   
        % estimating the energy     
        
        for l = delta:delta:1
            q1_minus = A * q0;
                
            %[yy, qq0, P1] = estimate_kf(A, B, C, u, q0, P0, Q, R, ...
            %    pow(k), k, eps);
        
            P1_minus = A * P0 * transpose(A) + Q;
    
            % estimation
            K1 = (P1_minus * transpose(C)) / ...
                (C * P1_minus * transpose(C) + R);
            q1 = q1_minus + K1 * (simpow - C * q1_minus);
            P1 = (eye(size(q0, 1)) + K1 * C) * P1_minus;
        
            y0_estimate = C * q1;
    
            % updating values for the next iteration
            q0 = q1;
            
            P0 = P1;
                
            y = [y; y0_estimate];
        
            q = [q; q0];
        end
                        
        k = k + 1;
    end
        
    i = i + 1;
end

%% plots

% plotting the path
subplot(2,2,[1 2]);
plot(pos(:, 1), pos(:, 2), 'Color', 'r', 'LineWidth', 1.2)

% plotting the energy and the estimated energy
subplot(2,2,3);
plot(linspace(0, k, size(pow, 1)), pow(:, 1))

subplot(2,2,4);
plot(linspace(0, k, size(y, 1)), y(:, 1))


%% model building function

function [A, C] = build_model(omega,r)

    m = 2*r+1;

    Aj = @(omega, j) [0 omega*j ; -omega*j 0];
    A  = zeros(m);
    A(1,1) = 0;
    C = [1];

    for i = 1:r
        A(2*i : 2*i + 1, 2*i : 2*i + 1) = Aj(omega, i); 
        C = [C 1 0];
    end
    
    C = 1/1 * C;

end
