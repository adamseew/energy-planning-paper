
% SIM6

% complete simulation of the algo (new physics, from Hector)



%% init


%%% initializing all the params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prompts with data for input

% asking data about the path

answer = inputdlg(...
    {'vehicle direction [deg]:','wind speed [m/s]:',...
     'wind direction [deg]:','start coordinate x [m]:','y [m]:',...
     'max power [W]:','min power [W]:','triggering point radius [m]'
    }, ...
    'path initialization',[1 40],...
    {'270','5','90','-100','220','60','30','10'});

if isempty(answer)
    strp = [270; 2; 345; -100; 220; 60; 30; 10];
else
    strp = str2double(answer);
end

clear answer;

% asking data about the algorithm

answer = inputdlg(...
    {'order r:','epsilon:','horizon N:'},...
    'model and algorithm initialization',[1 40],...
    {'3','1','60'});

if isempty(answer)
    strp2 = [3; 1; 60]; % default initial data
else
    strp2 = str2double(answer);
end

clear answer;


%%% gains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         .0005    .025 .001 
% kp, kvv, kd, ke1, ke2, ke3, ke4
strp3 = [5 5 .001 .0003 .02 .0006 .07];


%%% physics data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% storing data

vd = strp(1); % initial UAV direction
ws = strp(2); % wind speed
wd = strp(3); % wind direction
start_x = strp(4); % starting position x
start_y = strp(5); % y
max_pw = strp(6); % maximum reachable power
min_pw = strp(7); % minimum
trig_eps = [strp(8);strp(8)]; % radius of the triggering point (it's a 
                              % circular area)
r = strp2(1); % model size (the bigger the more precise)
eps = strp2(2); % optimization decrement step (to adapt control)
N = strp2(3); % optimization horizon (for MPC)

% dynamics control

m = 1; % mass of the aircraft [kg]
W = m*9.8; % weight force [N]

cl = 9.8/15^2; % coefficient to be determined experimentally
cth = 15/50; % same

kp = strp3(1); % positive gain constant to be also determined experimentally
kvv = strp3(2);
th_nominal = 50; % nominal vlalue of the throttle
th_delta = 0;
hd = 25; % desired altitude [m]

% guidance
kd = strp3(3);

% physics init
w =  ws*[cosd(wd);sind(wd)]; % wind vector (x, y axis velocity)

h = 20; % initial altitude
vv = 0; % initial vertical velocity
p = [start_x;start_y]; % position
theta = deg2rad(vd); % initial angle

sh = cth*(th_nominal+th_delta);
pdot = sh*[cos(theta);sin(theta)]+w; % initial velocity

wbx = dot(w,[cos(theta), sin(theta)]);
vs = cth*(th_nominal+th_delta)+wbx; % initial airspeed


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


%%% path plan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p1x_ = str2sym("115-sqrt(4900+c1)"); % start points first paths
p1y_ = str2sym("-146");
p2_ = str2sym("115");
p3x_ = str2sym("115-sqrt(5625)");
p3y_ = str2sym("-11");
p4_ = str2sym("115-2*sqrt(5625)");
psx_ = str2sym("-(2*75-2*sqrt(4900+c1))"); % shift
psy_ = str2sym("-(2*(75-sqrt(4900+c1))/10)");
ke1 = strp3(4); % gain first path
ke2 = strp3(5);
ke3 = strp3(6);
ke4 = strp3(7);


%%% params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beaware this works only with the default plan and has to be changed for a
% new one
min_c1 = -1000;
max_c1 = 0;
c1 = max_c1;
delta_params = 200;
n = 4; % used to measure the period T

psx = 0;
psy = 0;
p1x = double(subs(p1x_)); % to be redone every time c1 changes
p1y = double(subs(p1y_));
p2 = double(subs(p2_));
p3x = double(subs(p3x_));
p3y = double(subs(p3y_));
p4 = double(subs(p4_));
t1 = str2sym(strcat("(x+",num2str(p1x),")^2+(y+",num2str(p1y),...
    ")^2-4900-",num2str(c1)));
t2 = str2sym(strcat("x+",num2str(p2)));
t3 = str2sym(strcat("(x+",num2str(p3x),")^2+(y+",num2str(p3y),")^2-5625"));
t4 = str2sym(strcat("x+",num2str(p4)));
trg4 = [-p4;-p1y];
trg3 = [-p4;p3y];
trg2 = [-p2;p3y];
trg1 = [-p2;-p1y];


%%% loggers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_p = [start_x;start_y]; % position
log_y = []; % model output
log_q = []; % model state
log_period = []; % periods
log_pdangle = []; % vector field to pathectory angles
log_pow = []; % simulated power
log_dpd = []; % vector field result

log_vv = vv;
log_h = h;
log_theta = theta;
log_th_delta = th_delta;
log_vh = [];
log_pdot = pdot;

strp4 = [m W cl cth th_nominal th_delta hd h vv delta_T]; % saving 
                                                            % parameters


%% sim


k = 0;
i = 1;
time = 0;
v = 0;
syms x y;

last_trig = [175;161]; % the last triggering point
% next are just some simulation utilities
changed = 0; % for c1 parameter adaptation simulation
reached = 0; % reached the end of the simulation

                                                 
%%% physics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_C = 1; % also just simulation utility

while true 
    
    % getting the current pathectory
    if path_C == 4 % path 4
        
        path = t4;
        ke = ke4;
        dir = 1;
        trig = trg4;

        % updating the paths for the next iterations
        path_C = 0;
        p1x = double(subs(p1x_)); % to be redone every time c1 changes
        p3x = double(subs(p3x_));
        p4 = double(subs(p4_));
        
        psx__ = double(subs(psx_));
        psy__ = double(subs(psy_));
        psx = psx+psx__;
        psy = psy+psy__;
        t1 = str2sym(strcat("(x+",num2str(p1x+psx),")^2+(y+",...
            num2str(p1y+psy),")^2-4900-",num2str(c1)));
        t2 = str2sym(strcat("x+",num2str(p2+psx)));
        t3 = str2sym(strcat("(x+",num2str(p3x+psx),")^2+(y+",...
            num2str(p3y+psy),")^2-5625"));
        t4 = str2sym(strcat("x+",num2str(p4+psx)));
        trg4 = [-p4-psx;-p1y-psy];
        trg3 = [-p4-psx;p3y-psy];
        trg2 = [-p2-psx;p3y-psy];
        trg1 = [-p2-psx;-p1y-psy];

    elseif path_C == 3 % path 3
        
        path = t3;
        ke = ke3;
        dir = 1;
        trig = trg3;
    elseif path_C == 2 % path 2
        
        path = t2;
        ke = ke2;
        dir = -1;
        trig = trg2;
    else % path 1
        
        path = t1;
        ke = ke1;
        dir = 1;
        trig = trg1;
    end
    grad = gradient(path,[x y]); 
    hess = hessian(path,[x y]);
    grad = matlabFunction(grad); % getting function_handle (way faster then
                                 % sym
    hess = matlabFunction(hess);
    path = matlabFunction(path);
    
    updated_period = 0; % used to check if the period has been updated
        
    if and(and(mod(i-1,n) == 0,time > 1),updated_period == 0)
        
        period = time;
        [A C] = build_model(2*pi/period,r);
        Ad = A*delta_T+eye(2*r+1);  
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
       
        u_theta = gvf_control_2D(log_p(:,end).',... % guidance
            pdot,ke,kd,path,grad,hess,dir);
        
        th_delta = kp*(hd-h)-kvv*vv;
        wbx = dot(w,[cos(theta),sin(theta)]);
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
        theta = wrapToPi(theta); % normalizing between -pi and pi
        
        % log
        log_vv = [log_vv;vv];
        log_h = [log_h;h];
        log_theta = [log_theta;theta];
        log_th_delta = [log_th_delta;th_delta];
        log_p = [log_p p];
        log_pdot = [log_pdot pdot];
                 
        % produce a simulated energy value
        % first, let's calculate the angle between the two vectors
        %angle = acosd((pos_x*wind_x+pos_y*wind_y)/...
        %    (sqrt(pos_x^2+pos_y^2)*sqrt(wind_x^2+wind_y^2))...
        %             );
        % angle is nan? no windspeed then; so the contribution of the
        % wind is none
        %if isnan(angle)
        %    angle = 0;
        %end
            
        % if the angle is 0, lowest energy (tail wind). 180, heighest
        % (head wind)
        %log_pow = [log_pow;...
        %    min_pw+(max_pw-min_pw)*.5*abs(1-cosd(angle))*ws/vs]; 
                                             % winddir and windspeed effect
        
        %%% estimating the energy with KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %q1_minus = Ad*q0;
        %P1_minus = Ad*P0*Ad.'+Q; % fix the covariance, just sim
        %K1 = (P1_minus*C.')*(C*P1_minus*C.'+R)^-1;
        %q1 = q1_minus+K1*(log_pow(end)-C*q1_minus);
        %P1 = (eye(size(q0,1))+K1*C)*P1_minus;
        %y0_estimate = C*q1;
        %q0 = q1;
                
        %log_y = [log_y;y0_estimate];
        %log_q = [log_q q0];
                                
        k = k+1;
                
        % came over the last triggering point
        if all(abs(log_p(:,end)-trig(:)) <= trig_eps) 
            break;
        end
        
        if all(log_p(:,end)>last_trig) % checking if reached the last 
                                       % triggering point
            reached = 1;
            break;
        end
        
        % reached the battery striking point (just trying)
        if and(k*delta_T >= 200,changed == 0)
            plot(log_p(1,:),log_p(2,:),'Color','r','LineWidth',1.2)
            
            % forcing params (just testing)
            %max_c1 = -1000;
            c1 = max_c1;
            
            p1x = double(subs(p1x_)); % to be redone every time c1 changes
            p3x = double(subs(p3x_));
            p4 = double(subs(p4_));
            
            psx = psx-psx__;
            psy = psy-psy__;
            psx__ = double(subs(psx_));
            psy__ = double(subs(psy_));
            psx = psx+psx__;
            psy = psy+psy__;
            t1 = str2sym(strcat("(x+",num2str(p1x+psx),")^2+(y+",...
                num2str(p1y+psy),")^2-4900-",num2str(c1)));
            t2 = str2sym(strcat("x+",num2str(p2+psx)));
            t3 = str2sym(strcat("(x+",num2str(p3x+psx),")^2+(y+",...
                num2str(p3y+psy),")^2-5625"));
            t4 = str2sym(strcat("x+",num2str(p4+psx)));
            trg4 = [-p4-psx;-p1y-psy];
            trg3 = [-p4-psx;p3y-psy];
            trg2 = [-p2-psx;p3y-psy];
            trg1 = [-p2-psx;-p1y-psy];
            
            if path_C == 4 % updating paths with new data
                path = t4;
                trig = trg4;
            elseif path_C == 3
                path = t3;
                trig = trg3;
            elseif path_C == 2
                path = t2;
                trig = trg2;
            else
                path = t1;
                trig = trg1;
            end
            grad = gradient(path,[x y]); 
            hess = hessian(path,[x y]);
            grad = matlabFunction(grad);
            hess = matlabFunction(hess);
            path = matlabFunction(path);
            changed = 1;
        end
    end
    
    if reached
        break;
    end
    
    i = i+1;
    path_C = path_C+1;
    
end  



%% plots


figure(1);
time_v = linspace(0,k*delta_T,length(log_pow));
% plotting the path
subplot(4,2,[1 2]);
plot(log_p(1,:),log_p(2,:),'Color','r')
title('trajectory');
% plotting the energy and the estimated energy
subplot(4,2,5);
plot(time_v,log_pow(:,1))
title('energy sensor');
subplot(4,2,6);
plot(time_v,log_y)
title('estimated energy');
subplot(4,2,3);
plot(time_v,log_h);
title('altitute');
subplot(4,2,4);
plot(time_v,log_th_delta+th_nominal)
title('throttle');

figure(2); % to be fixed if size is different from 3
subplot(7,1,1);
plot(time_v,log_q(1,1:end))
title('model coef alpha 0');
subplot(7,1,2);
plot(time_v,log_q(2,1:end))
title('alpha 1');
subplot(7,1,3);
plot(time_v,log_q(3,1:end))
title('beta 1');
subplot(7,1,4);
plot(time_v,log_q(4,1:end))
title('alpha 2');
subplot(7,1,5);
plot(time_v,log_q(5,1:end))
title('beta 2');
subplot(7,1,6);
plot(time_v,log_q(6,1:end))
title('alpha 3');
subplot(7,1,7);
plot(time_v,log_q(7,1:end))
title('beta 3');



%% save


% asking data about the algorithm

answer = inputdlg(...
    {'simulation data label:',[1 40],...
    ''});

if isempty(answer)
    strp5 = 'NAME'; % default label
else
    strp5 = answer;
end

clear answer;

save(strcat(strp5,'.mat'));

csvwrite(strcat('position_simulation',strp5,'.csv'),[log_p(1,2:end).'...
    log_p(2,2:end).' log_pdangle log_dpd(:,1) log_dpd(:,2)]);
csvwrite(strcat('energy_simulation',strp5,'.csv'),[time_v' ...
    log_pow log_y log_q(1,:).' log_q(2,:).' log_q(3,:).' log_q(4,:).'...
    log_q(5,:).' log_q(6,:).' log_q(7,:).']);
csvwrite(strcat('data_simulation',strp5,'.csv'),[strp strp2 strp3 strp4]);
csvwrite(strcat('perioddata_simulation',strp5,'.csv'),log_period);



%% functions


function [u_theta] = gvf_control_2D(p,dot_p,ke,kd,path,grad,hess,dir)
%GVF_CONTROL_2D

    if (nargin(path) > 1) % function handle has two arguments (circle)
        e = path(p(:,1),p(:,2));    
        n = grad(p(:,1),p(:,2));
    else % one argument (line)
        e = path(p(:,1));    
        n = grad();
    end
    
    H = hess();
    E = [0 -1;1 0];
    tau = dir*E*n;

    dot_pd = tau-ke*e*n; % (7)
    ddot_pd = (E-ke*e)*H*dot_p-ke*n'*dot_p*n; % (10)
    ddot_pdhat = -E*(dot_pd*dot_pd')*E*ddot_pd/norm(dot_pd)^3; % (9)

    dot_Xid = ddot_pdhat'*E*dot_pd/norm(dot_pd); % (13)

    u_theta = dot_Xid+kd*dot_p'*E*dot_pd/(norm(dot_p)*norm(dot_pd)); % (16)

end

function [A, C] = build_model(omega,r)
%BUILD_MODEL Creates the simplified model (local implementation from
%build_model)
    m = 2*r+1;

    Aj = @(omega, j) [0 omega*j;-omega*j 0];
    A  = zeros(m);
    A(1,1) = 0;
    C = [1];

    for i = 1:r
        A(2*i:2*i+1,2*i:2*i+1) = Aj(omega,i); 
        C = [C 1 0];
    end

end

