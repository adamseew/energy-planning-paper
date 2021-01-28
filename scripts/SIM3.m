
% SIM3


%% init

%%% initializing all the params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prompts with data for input

% asking data about the path

answer = inputdlg(...
    {'horizontal speed [m/s]:','wind speed [m/s]:',...
     'vehicle direction [angles]:','wind direction [angles]:',...
     'start coordinate x [m]:','start coordinate y [m]:',...
     'max power [W]:','min power [W]:','triggering point radius [m]',...
     'altitude [m]'...
    }, ...
    'path initialization',[1 40],...
    {'20','5','270','90','-100','220','60','30','10', '25'}); 

if isempty(answer)
    strp = [20; 5; 270; 90; -100; 220; 60; 30; 10; 25]; % default
                                                        % initial data
else
    strp = str2double(answer);
end

clear answer;

% asking data about the algorithm

answer = inputdlg(...
    {'order r:','epsilon:','horizon N:'},...
    'model and algorithm initialization',[1 40],...
    {'3','1','10'});

if isempty(answer)
    strp2 = [3; 1; 10]; % default initial data
else
    strp2 = str2double(answer);
end

clear answer;

% stroing data

delta_T = 1e-2;

m = 1; % mass of the aircraft [kg]
W = m*9.8; % weight force [N]

cl = 9.8/15^2; % coefficient to be determined experimentally
cth = 15/50; % same

kp = 5; % positive gain constant to be also determined experimentally
kvv = 5;
th_nominal = 50; % nominal vlalue of the throttle
th_delta = 0;
hd = 25; % desired altitude [m]

w = [strp(2)*cosd(strp(4));strp(2)*sind(strp(4))]; % wind vector 
                                                   % (x, y axis velocity)

vv = 0; % vertical velocity
sh = strp(1); % horizontal speed
h = strp(10); % altitude
p = [strp(5); strp(6)]; % position
theta = deg2rad(strp(3)); % initial angle
vh = sh * [cos(theta);sin(theta)]; % initial velocity
wbx = dot(w,[cos(theta), sin(theta)]);
vs = cth*(th_nominal + th_delta) + wbx; % initial airspeed

maxpw = strp(7); % minimum power of the drone
minpw = strp(8); % maximum
trigeps = strp(9); % radius of the triggering points

r = strp2(1); % number of coefficients in the modle (more means more 
              % accuracy)
eps = strp2(2); % MPC related
N = strp2(3); % MPC horizon



%%% building the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

period = 1;
[A C] = build_model(2*pi/1, r);
syms q [2*r+1 1];
% discretizing 
Ad = A*delta_T+eye(2*r+1);
% just a very approximate guesses!
q0 = ones(size(q,1),1); % initial state guess
P0 = ones(size(q0,1)); % covariance guess
% process noise and sensor noise
Q = ones(size(q,1));
R = 1;


% asking for the plan specifization (cancel to get the default one

disp('plan specification');
[bn folder] = uigetfile('.pln');
if bn == 0
    fid = fopen('../data/simulation3/plan_specification.pln','rt');
else
    fid = fopen(fullfile([folder bn]),'rt');
end



%%% parsing plan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = []; % trajectories

trigs = []; % triggering points

i = -1;
shift = [];
while true
    thisline = fgetl(fid);
    if ~ischar(thisline) 
        break; 
    elseif isempty(thisline)
        continue;
    end  % end of file
    
    if i == -1 % first line contains the shift
        thisline = split(thisline,',');
        shift = [str2double(thisline(1)); str2double(thisline(2))];
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



%%% iterating trajectories to get n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = [];
p_fix = [0; 0];
n = 0;

% n is used to measure the period T

for traj = transpose(path) % per each trajectory
        
    traj = split(traj,";");
    traj = str2sym(traj(3));
    x = p_fix(1);
    y = p_fix(1);
    e_j = double(subs(traj));
    x = p_fix(1)-shift(1);
    y = p_fix(2)-shift(2);
    
    if ismember(double(subs(traj)),e)
        break;
    else
        e = [e e_j];
    end
    
    n = n+1;
end

% fprintf('n is %d\n', n);



%%% loggers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_vv = [];
log_h = [];
log_theta = [];
log_sh = log_vv;
log_th_delta = log_vv;
log_p = [];
log_vh = [];

log_vv = [log_vv; vv];
log_h = [log_h; h];
log_theta = [log_theta; theta];
log_sh = [log_sh; sh];
log_th_delta = [log_th_delta; th_delta];
log_p = [log_p; p.'];
log_vh = [log_vh vh];

log_pdangle = []; % vector field
log_y = []; % model output
log_q = []; % model state
log_pow = []; % simulated power
log_period = [];



%% sim

addpath(genpath('position')); % needed for build_gdn2

k = 0;
i = 1;

time = 0;

for traj = transpose(path)
        
    traj = split(traj, ";"); % get the i-th trajectory from the plan
    updated_period = 0; % used to check if the period has been updated
    if contains(traj(2), '1') % rotation (E)
        E = [0 1; -1 0];
    else
        E = [0 -1; 1 0];
    end
        
    ke = str2double(traj(1)); % gain (speed of convergence to the 
                              % trajectory)
    
    if and(and(mod(i-1,n) == 0,time > 1),updated_period == 0)
        
        period = time;
        log_period = [log_period; k*delta_T period];
        [A C] = build_model(2*pi/period,r); % re-initialize model
        Ad = A*delta+eye(2*r+1); % discretize
        time = 0; % restart measuring time to check whenever it changed at
                  % the next period
        updated_period = 1;
    end
    
    while true
        
        time = time+delta_T;        
            
        if and(period < time,time > 1)
                   
            period = time;
            log_period = [log_period; k*delta_T period];
            [A C] = build_model(2*pi/period, r);
            Ad = A*delta_T+eye(2*r+1);
        end
            
        %%% vector field update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        [dpd, pdangle] = build_gdn2(E,ke,str2sym(traj(3)),p.'); 
                        
        log_pdangle = [log_pdangle; pdangle];
        
        % is the control u pdangle? NOPE
        u_theta=deg2rad(pdangle); % radians 
        
        % Vertical dynamics
        th_delta = kp*(hd-h)-kvv*vv;
        wbx = dot(w,[cos(theta) sin(theta)]);
        vs = cth*(th_nominal+th_delta)+wbx;
        L = cl*vs*vs;
        av = 1/m*(L-W);

        vv = vv+av*delta_T;
        h = h+vv*delta_T;

        % Horizontal unicycle model
        sh = cth*(th_nominal+th_delta);
        pdot = sh*[cos(theta);sin(theta)]+w;

        p = p+pdot*delta_T;
        theta = theta+u_theta*delta_T;
        
        log_vv = [log_vv; vv];
        log_h = [log_h; h];
        log_theta = [log_theta; theta];
        log_sh = [log_sh; sh];
        log_th_delta = [log_th_delta; th_delta];
        log_p = [log_p; p.'];
        log_vh = [log_vh vh];
          
        
        
        %%% energy simulation data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % gives a simulated energy value
        % throttle goes from nominal +- nominal (e.g., if the nominal
        % throttle is 50, it goes from 0 to 100)
        % then we just scale it with the energy
        log_pow = [log_pow; minpw+...
                   (th_nominal+th_delta/(2*th_nominal))*(maxpw-minpw)
                  ]; % energy signal
        
        
        
        %%% KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        % estimating the energy     
        
        q1_minus = Ad*q0;
        
        P1_minus = Ad*P0*Ad.'+Q; % fix the covariance, just sim
        % estimation
        K1 = (P1_minus*C.')*(C*P1_minus*C.'+R)^-1;
        q1 = q1_minus+K1*(log_pow(end)-C*q1_minus);
        P1 = (eye(size(q0,1))+K1*C)*P1_minus;
        y0_estimate = C*q1;
        % updating values for the next iteration
        q0 = q1;
        log_y = [log_y; y0_estimate];
        log_q = [log_q q0];
        
        k = k+1;
        
        % reached the triggering point, going to stage i+1
        if all(abs(p - trigs(i,:)) <= trigeps) 
            break;
        end
    end
        
    i = i+1;
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


%% save

csvwrite('position_simulationNAME.csv', [pos(:, 1) pos(:, 2) hlist pdanglelist]);

csvwrite('energy_simulationNAME.csv', [linspace(0, size(pow, 1)*delta, size(pow, 1))', ...
    pow(:, 1), y(:, 1), q(1,1:end)', q(2,1:end)', q(3,1:end)',     ...
    q(4,1:end)', q(5,1:end)', q(6,1:end)', q(7,1:end)']);

csvwrite('trajdata_simulationNAME.csv', strp);

csvwrite('algdata_simulationNAME.csv', strp2);

csvwrite('perioddata_simulationNAME.csv', periodlist);


%% model

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

