
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

% asking data about the algorithm

answer = inputdlg(...
    {'order r:','epsilon:','horizon N:'},...
    'model and algorithm initialization',[1 40],...
    {'3','1','10'}); % asking for initial data

if isempty(answer)
    strp2 = [3; 1; 10]; % default initial data
else
    strp2 = str2double(answer);
end

clear answer;

% stroing data

r = strp2(1);
eps = strp2(2);
N = strp2(3);

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

n = i;

fclose(fid);


%% building the model

period = 1;

delta = .01;
[A C] = build_model(2*pi/1, r);
syms q [2*r+1 1];

% discretizing 

Ad = A*delta+eye(2*r+1);
%Ad = exp(A*delta);

%% model guesses

q0 = ones(size(q,1),1);

P0 = ones(size(q,1));
    
% process noise and sensor noise

Q = ones(size(q,1));
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
pown = [];    

nowpos = [startx starty];
pdanglelist = [];
    
% contribution of the wind. This doesn't change as winspeed and
% direction is constant
windx = delta * windspeed * cosd(winddir);  % using delta to make the path more smooth
windy = delta * windspeed * sind(winddir);
    
y = []; % model output
q = []; % model state

addpath(genpath('position')); % needed for build_gdn2

i = 1;
time = 0;

v = 0;

periodlist = [];

for traj = transpose(path)
        
    traj = split(traj, ";");
            
    if contains(traj(2), '1')
        E = [0 1; -1 0];
    else
        E = [0 -1; 1 0];
    end
        
    ke = str2double(traj(1));
    
    while true
        
        time = time + delta;
        
        if period < time
        
            period = time;
            
            periodlist = [periodlist; k period];
        
            [A C] = build_model(2*pi/period, r);
            Ad = A*delta+eye(2*r+1);
        
            % re-initialization of the KF
            P0 = ones(size(q,1));
            
        elseif mod(i, n) == 0
        
            period = time;
            
            periodlist = [periodlist; k period];
        
            [A C] = build_model(2*pi/period, r);
            Ad = A*delta+eye(2*r+1);
            
            % re-initialization of the KF
            P0 = ones(size(q,1));
        
            time = 1;
        
        end

        pos = [pos; nowpos];
            
        % vector field
        [dpd, pdangle] = build_gdn2(E, ke, str2sym(traj(3)), nowpos); 
                        
        pdanglelist = [pdanglelist; pdangle];
        
        % dpd is the ideal direction... We have wind though
        % now if the time is sampled every second, the below expression
        % actually indicates the offset from the original location
            
        posx = delta * vehspeed * cosd(pdangle);
        posy = delta * vehspeed * sind(pdangle);
        
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
        
        %if mod(k, 10) == 0
        %    v = randn; % uncomment to add some noise
        %end
        %simpow = simpow+v; % energy signal with normal noise
        
        pow = [pow; simpow]; % energy signal
        
                   
        % estimating the energy     
        q1_minus = Ad * q0;
        
        P1_minus = Ad * P0 * transpose(Ad) + Q;
        
        % estimation
        K1 = (P1_minus * transpose(C)) / ...
            (C * P1_minus * transpose(C) + R);
        q1 = q1_minus + K1 * (pow(end) - C * q1_minus);
        P1 = (eye(size(q0, 1)) + K1 * C) * P1_minus;
        
        y0_estimate = C * q1;
    
            % updating values for the next iteration
        q0 = q1;
            
        P0 = P1;
                
        y = [y; y0_estimate];
        
        q = [q q0];
                                
        k = k + 1;
    end
        
    i = i + 1;
end


%% plots

figure(1);

% plotting the path
subplot(2,2,[1 2]);
plot(pos(:, 1), pos(:, 2), 'Color', 'r', 'LineWidth', 1.2)

% plotting the energy and the estimated energy
subplot(2,2,3);
plot(linspace(0, k, size(pow, 1)), pow(:, 1))

subplot(2,2,4);
plot(linspace(0, k, size(y, 1)), y(:, 1))

figure(2);
subplot(7,1,1);
plot(q(1,1:end))

subplot(7,1,2);
plot(q(2,1:end))

subplot(7,1,3);
plot(q(3,1:end))

subplot(7,1,4);
plot(q(4,1:end))

subplot(7,1,5);
plot(q(5,1:end))

subplot(7,1,6);
plot(q(6,1:end))

subplot(7,1,7);
plot(q(7,1:end))


%% saving data

csvwrite('position_simulation3A.csv', [pos(:, 1) pos(:, 2) pdanglelist]);

csvwrite('energy_simulation3A.csv', [linspace(0, k, size(pow, 1))', ...
    pow(:, 1), y(:, 1), q(1,1:end)', q(2,1:end)', q(3,1:end)',     ...
    q(4,1:end)', q(5,1:end)', q(6,1:end)', q(7,1:end)']);

csvwrite('trajdata_simulation3A.csv', strp);

csvwrite('algdata_simulation3A.csv', strp2);

csvwrite('perioddata_simulation3A.csv', periodlist);


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

end
