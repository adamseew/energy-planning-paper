
% SIM5

% complete simulation of the algo (old physics, not from Hector)



%% init


%%% initializing all the params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prompts with data for input

% asking data about the path

answer = inputdlg(...
    {'vehicle speed [m/s]:','wind speed [m/s]:',...
     'vehicle direction [angles]:','wind direction [angles]:',...
     'start coordinate x [m]:','start coordinate y [m]:',...
     'max power [W]:','min power [W]:','triggering point radius [m]'
    }, ...
    'path initialization',[1 40],...
    {'20','5','270','90','-100','220','60','30','10'});

if isempty(answer)
    strp = [20; 5; 270; 90; -100; 220; 60; 30; 10];
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

% storing data

vs = strp(1); % desired vecrtical speed
ws = strp(2); % wind speed
vd = strp(3); % initial UAV direction 
wd = strp(4); % wind direction
start_x = strp(5); % starting position x
start_y = strp(6); % y
max_pw = strp(7); % maximum reachable power
min_pw = strp(8); % minimum
trig_eps = strp(9); % radius of the triggering point (it's a circular area)

r = strp2(1); % model size (the bigger the more precise)
eps = strp2(2); % optimization decrement step (to adapt control)
N = strp2(3); % optimization horizon (for MPC)



%%% building the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_T =  1e-2;
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


% asking for the plan specifization (cancel to get the default one)

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
        shift = [str2double(thisline(1));str2double(thisline(2))];
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



%%% iterating trajectories to get n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (to measure the period)

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



%%% loggers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_p = [start_x;start_y];
log_y = [];
log_q = [];
log_period = [];
log_pdangle = [];
log_pow = [];
log_dpd = [];
 


%% sim
    
addpath(genpath('position')); % needed for build_gdn2

k = 0;
i = 1;
time = 0;
v = 0;

% contribution of the wind. This doesn't change as winspeed and
% direction are constant
wind_x = delta_T*ws*cosd(wd); 
wind_y = delta_T*ws*sind(wd);
y = []; % model output
q = []; % model state

for traj = transpose(path)
        
    traj = split(traj,";");
    updated_period = 0;
            
    if contains(traj(2),'1')
        E = [0 1;-1 0];
    else
        E = [0 -1;1 0];
    end
        
    ke = str2double(traj(1));
    
    if and(and(mod(i-1,n) == 0,time > 1),updated_period == 0)

        period = time;
        [A C] = build_model(2*pi/period,r);
        Ad = A*delta+eye(2*r+1);  
        time = 0;
        updated_period = 1;
    end
    
    while true
        
        time = time+delta_T;
        
        if and(period < time,time > 1)
                   
            period = time;
            log_period = [log_period; k period];
            [A C] = build_model(2*pi/period,r);
            Ad = A*delta_T+eye(2*r+1);
        end
        
        %%% vector field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [dpd, pdangle] = build_gdn2(E,ke,str2sym(traj(3)),log_p(:,end).');
        log_pdangle = [log_pdangle;pdangle];
        log_dpd = [log_dpd dpd];
        
        % dpd is the ideal direction... We have wind though
        % now if the time is sampled every second, the below expression
        % actually indicates the offset from the original location
        pos_x = delta_T*vs*cosd(pdangle);
        pos_y = delta_T*vs*sind(pdangle);
        
        % new position
        log_p = [log_p ...
            [log_p(1,end)+wind_x+pos_x;log_p(2,end)+wind_y+pos_y]]; 
                 
        % produce a simulated energy value
        % first, let's calculate the angle between the two vectors
        angle = acosd((pos_x*wind_x+pos_y*wind_y)/...
            (sqrt(pos_x^2+pos_y^2)*sqrt(wind_x^2+wind_y^2))...
                     );
        % angle is nan? no windspeed then; so the contribution of the
        % wind is none
        if isnan(angle)
            angle = 0;
        end
            
        % if the angle is 0, lowest energy (tail wind). 180, heighest
        % (head wind)
        log_pow = [log_pow;...
            min_pw+(max_pw-min_pw)*.5*abs(1-cosd(angle))*ws/vs]; 
                                             % winddir and windspeed effect
        
        %%% estimating the energy with KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        q1_minus = Ad*q0;
        P1_minus = Ad*P0*Ad.'+Q; % fix the covariance, just sim
        K1 = (P1_minus*C.')*(C*P1_minus*C.'+R)^-1;
        q1 = q1_minus+K1*(log_pow(end)-C*q1_minus);
        P1 = (eye(size(q0,1))+K1*C)*P1_minus;
        y0_estimate = C*q1;
        q0 = q1;
                
        log_y = [log_y;y0_estimate];
        log_q = [log_q q0];
                                
        k = k+1;
        
        % reached the triggering point, going to stage i+1
        if all(abs(log_p(end)-trigs(i,:)) <= trig_eps) 
            break;
        end
    end
    
    i = i+1;
end



%% plots

figure(1);

time_v = linspace(0,k,size(log_p,1));
% plotting the path
subplot(2,2,[1 2]);
plot(log_p(:,1),log_p(:,2),'Color','r','LineWidth',1.2)
% plotting the energy and the estimated energy
subplot(2,2,3);
plot(time_v,log_pow(:,1))
subplot(2,2,4);
plot(time_v,log_y)

figure(2); % to be fixed if size is different from r=3
subplot(7,1,1);
plot(time_v,log_q(1,1:end))
subplot(7,1,2);
plot(time_v,log_q(2,1:end))
subplot(7,1,3);
plot(time_v,log_q(3,1:end))
subplot(7,1,4);
plot(time_v,log_q(4,1:end))
subplot(7,1,5);
plot(time_v,log_q(5,1:end))
subplot(7,1,6);
plot(time_v,log_q(6,1:end))
subplot(7,1,7);
plot(time_v,log_q(7,1:end))



%% save

csvwrite('position_simulationNAME.csv',[log_p(:,1).' log_p(:,2).' ...
    log_pdangle log_dpd(:,1).' log_dpd(:,2).']);
csvwrite('energy_simulationNAME.csv',[time_v', ...
    log_pow(:,1),log_y(:,1),log_q(1,1:end).',log_q(2,1:end).',...
    log_q(3,1:end).',log_q(4,1:end).',log_q(5,1:end).',log_q(6,1:end).',...
    log_q(7,1:end).']);
csvwrite('trajdata_simulationNAME.csv',strp);
csvwrite('algdata_simulationNAME.csv',strp2);
csvwrite('perioddata_simulationNAME.csv',log_period);



%% functions

function [A, C] = build_model(omega,r)
%BUILD_MODEL Creates the simplified model (local implementation from
%build_model)
    m = 2*r+1;

    Aj = @(omega, j) [0 omega*j ; -omega*j 0];
    A  = zeros(m);
    A(1,1) = 0;
    C = [1];

    for i = 1:r
        A(2*i : 2*i+1, 2*i : 2*i+1) = Aj(omega,i); 
        C = [C 1 0];
    end

end


