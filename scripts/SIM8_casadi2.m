
% SIM8

% simulation of the algorithm, path I, dynamic, revised MPC


%% casadi opt control

addpath('../../../../../');
% UPDATE THE PATH WITH CASADI SOURCE !


%% cleaning

clear
clc


%% init


%%% initializing all the params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prompts with data for input

% asking data about the path
% uncomment to manualy change conditions
%answer = inputdlg(...
%    {'vehicle direction [deg]:','wind speed [m/s]:',...
%     'wind direction [deg]:','start coordinate x [m]:','y [m]:',...
%     'max power [W]:','min power [W]:','triggering point radius [m]'
%    }, ...
%    'path initialization',[1 40],...
%    {'270','5','0','-100','220','36','16','20'});

%if isempty(answer)
strp = [270; 5; 0; -100; 220; 36; 16; 20];
%else
%    strp = str2double(answer);
%end

%clear answer;


%%% physics data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gains
%       kp, kvv, kd,  ke1, ke2, ke3, ke4
strp3 = [5   5  .001 .006  .1  .006  .1];
%         

vd = strp(1); % initial UAV direction
ws = strp(2); % wind speed
wd = strp(3); % wind direction
start_x = strp(4); % starting position x
start_y = strp(5); % y
max_pw = strp(6); % maximum reachable power
min_pw = strp(7); % minimum
trig_eps = [strp(8);strp(8)]; % radius of the triggering point (it's a 
                              % circular area)

strp2 = 3;
r = strp2; % model size (the bigger the more precise)

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

% battery
int_v = 14.8; % internal battery voltage [V]
qc = 3.2; % constance nominal capacity [Ah]
res = 0.0012; % resistance [ohm]


%%% building the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_T =  1e-2; % 1/100 of a second
period = 1;
[A B C] = build_model(2*pi/1, r);
syms q [2*r+1 1];
% discretizing 
Ad = A*delta_T+eye(2*r+1);
% just a very approximate guesses!
q0 = ones(size(q,1),1); % initial state guess
P0 = ones(size(q0,1)); % covariance guess
% process noise and sensor noise
Q = ones(size(q,1));
R = 1;
N = 6; % MPC horizon (6 seconds)

% soc
kb = .00183; % battery coefficient
b = @(y) -kb*(int_v-sqrt(int_v^2-4*res*y))/(2*res*qc); % battery dynamics
soc = linspace(1,.9,93);
soc = [soc linspace(.84,.45,177)];
soc = [soc linspace(.36,.2,730)]; % with sudden battery drops


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
% to be changed for a different plan from the one in the paper

min_c1 = -1000;
max_c1 = 0;
c1 = max_c1; % trajectory parameter
c1_old = c1;

min_c2 = 2;
max_c2 = 10;
c2 = max_c2; % computational parameter
c2_old = c2;

% scaling factors
max_t = 705.96; % time with max_c1
min_t = 270.15;
g_max_c2 = 9.75462; % g(max_c2)
g_min_c2 = 4.75376; % g(min_c2)

nu1 = -(max_t-min_t)/min_c1; % traj param
tau1 = -min_c1*(min_t-max_t)/min_c1+min_t;
nu2 = (g_max_c2-g_min_c2)/(max_c2-min_c2); % computational param
tau2 = min_c2*(g_min_c2-g_max_c2)/(max_c2-min_c2)+g_min_c2;


% control
est_u = @(c1,c2) diag([nu1 nu2])*([c1;c2])+[tau1;tau2];
est_u_old = [0;0];
u = @(est_u,est_u_old) est_u-est_u_old;

% delta for param controller
deltas = [500;2];

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
log_pow = []; % simulated power

log_vv = vv;
log_h = h;
log_theta = theta;
log_th_delta = th_delta;
log_pdot = pdot;

log_c1 = [];
log_c2 = [];

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
reached = 0; % reached the end of the simulation
path_C = 1; % also just simulation utility
get_optcotrol = 1; % inhibits the controller whenever the limit is reached

                                                 
%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while true 
    
    % getting the current pathectory
    if path_C == 4 % path 4
        
        path = t4;
        ke = ke4;
        dir = 1;
        trig = trg4;

        % updating the paths for the next iterations (to be shifted from
        % primitive paths
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
                                 % sym)
    hess = matlabFunction(hess);
    path = matlabFunction(path);
    
    updated_period = 0; % used to check if the period has been updated
        
    if and(and(mod(i-1,n) == 0,time > 1),updated_period == 0)
        
        period = time;
        [A B C] = build_model(2*pi/period,r);
        Ad = A*delta_T+eye(2*r+1);  
        time = 0;
        updated_period = 1;
    end
    
    while true
        
        time = time+delta_T;
        
        if and(period < time,time > 1)
                   
            period = time;
            log_period = [log_period;k period];
            [A B C] = build_model(2*pi/period,r);
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
                 
        % produce a simulated energy value from throttle
        log_pow = [log_pow;...
            (max_pw-min_pw)*(th_delta+th_nominal)/(2*th_nominal)+min_pw];
        
        
        %%% params controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if and(get_optcotrol*k*delta_T >= 60,... % don't optimize for 1 min
               get_optcotrol*k*delta_T >= 2*period) % for 2 periods anyways
            
            c2_chain = []; % chain of predicted opt controls from MPC
            
            if mod(k,1/delta_T) == 0 % MPC runs every second
                
                c2_chain = mpc(min_c1,max_c1,min_c2,max_c2,c1,c2,N,eu,...
                    soc(round(k*delta_T)+1),b,qc,int_v,q0,Ad,B,C,u,...
                    est_u,k,delta_T);
                
                old_c2 = c2;
                c2 = c2_chain(2); % getting c2
                
            end
                
            % estimating remaining time from c1
            sim_t = est_u(c1,c2);
            sim_t = sim_t(1); % remaining time if the param was changed
                              % at the beginning
            rem_t = (sim_t/max_t)*(max_t-k*delta_T); % at the current 
                                                     % instant
            
            b0 = soc(round(k*delta_T)+1); % battery state
            qq0 = q0; % state at this time instant                                         
            bat_t = 0; % battery time                                         
                                                     
            eeu = eu; % control estimate (from model)
            eeeu =  est_u(c1,c2); 
            
            c2_mpc_count = 2;
            
            while b0 > 0
                
                if and(~isempty(c2_chain),...
                       length(c2_chain) >= c2_mpc_count)
                        eeeu = est_u(c1,c2_chain(c2_mpc_count));
                        c2_mpc_count = c2_mpc_count+1;
                end
                
                for jj = delta_T:delta_T:1 % 1 sec evolution with euler
                    qq1 = Ad*qq0+B*u(eeeu,eeu);
                    yy1 = C*qq1;
                    b0 = b0+delta_T*b(yy1);
                    qq0 = qq1;
                    eeu = eeeu;
                    bat_t = bat_t+delta_T;
            
                    if bat_t > rem_t
                        break; % all good
                    end
                end
                    
                if bat_t > rem_t
                    break; % all good (it reached this point from above)
                end
            end
                
            if rem_t > bat_t
                if c1-deltas(1) >= min_c1
                    c1 = c1-deltas(1); % reducing path!
                end
            end
            
            if and(c1 <= min_c1,c2 <= min_c2) % lowest possible control 
                                              % reached (this changes in
                                              % path II where you start 
                                              % from lowest and go to the
                                              % highest)
                get_optcotrol = 0;
            end
        end
        
        %%% estimating the energy with KF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        eu = est_u(c1,c2); % u estimate
        log_pow(end) = log_pow(end)+eu(2); % eu(2) simulates the computa-
                                           % tional power
        q1_minus = Ad*q0+B*u(eu,est_u_old);
        est_u_old = eu;
        P1_minus = Ad*P0*Ad.'+Q;
        K1 = (P1_minus*C.')*(C*P1_minus*C.'+R)^-1;
        q1 = q1_minus+K1*(log_pow(end)-C*q1_minus);
        P1 = (eye(size(q0,1))+K1*C)*P1_minus;
        y0_estimate = C*q1;
        q0 = q1;
                
        log_y = [log_y;y0_estimate];
        log_q = [log_q q0];
        log_c1 = [log_c1;c1];
        log_c2 = [log_c2;c2];
                                
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
        
        % traj param changed, we need to change the params
        if c1_old ~= c1
            
            p1x = double(subs(p1x_));
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
        end
        
        c1_old = c1; % saving params for next iteration
        c2_old = c2;
        
    end
    
    if reached
        break;
    end
    
    i = i+1;
    path_C = path_C+1;
    
end  

time_v = linspace(0,k*delta_T,length(log_pow)); % time vector
log_b = soc(1:round(k*delta_T)+1);
log_b = [linspace(0,k*delta_T,length(log_b)).' log_b.' log_b.'*qc*int_v];


%% plots


figure; % plotting the path
plot(log_p(1,:),log_p(2,:),'Color','r')
title('model traj');
xlabel('x (m)');
ylabel('y (m)');

figure; % plotting the energy and the estimated energy
t = tiledlayout(2,2);
nexttile,plot(time_v,log_pow(:,1))
title('energy sensor');
ylabel('power (W)');
nexttile,plot(time_v,log_y)
title('estimated energy');
ylabel('power (W)');
nexttile,plot(time_v,log_h(2:end));
title('altitute');
ylabel('z (m)');
nexttile,plot(time_v,log_th_delta(2:end)+th_nominal)
title('throttle');
ylabel('value');
title(t,'model energy')
xlabel(t,'time (sec)')

figure;
plot(time_v,log_y)
xlabel('time (sec)')
ylabel('power (W)')
title('algorithms')
hold on
plot(log_b(:,1),log_b(:,3))
axes('Position',[.6 .6 .3 .3])
box on
plot(time_v,log_c1)
ylabel('c1')
hold on
yyaxis right
plot(time_v,log_c2)
ylabel('c2')

figure; % to be fixed if size is different from 3
t = tiledlayout(7,1);
nexttile,plot(time_v,log_q(1,1:end))
title('alpha 0');
nexttile,plot(time_v,log_q(2,1:end))
title('alpha 1');
nexttile,plot(time_v,log_q(3,1:end))
title('beta 1');
nexttile,plot(time_v,log_q(4,1:end))
title('alpha 2');
nexttile,plot(time_v,log_q(5,1:end))
title('beta 2');
nexttile,plot(time_v,log_q(6,1:end))
title('alpha 3');
nexttile,plot(time_v,log_q(7,1:end))
title('beta 3');
title(t,'model coefs')
xlabel(t,'time (sec)')
ylabel(t,'value')


%% save


% asking data about the algorithm

answer = inputdlg('simulation data label:','save',[1 40],{'NAME'});

if isempty(answer)
    clear answer;
    return; % pressed cancel, do not save
else
    strp5 = answer;
end

clear answer;

strp5 = string(strp5);
save(strcat(strp5,'.mat'));

csvwrite(strcat('position_simulation',strp5,'.csv'),[log_p(1,:).'...
    log_p(2,:).' log_h log_theta log_vv log_th_delta ...
    log_pdot(1,:).' log_pdot(2,:).']);
csvwrite(strcat('energy_simulation',strp5,'.csv'),[time_v.' ...
    log_pow log_y log_q(1,:).' log_q(2,:).' log_q(3,:).' log_q(4,:).'...
    log_q(5,:).' log_q(6,:).' log_q(7,:).']);
csvwrite(strcat('perioddata_simulation',strp5,'.csv'),log_period);
csvwrite(strcat('data_simulation',strp5,'.csv'),[strp.' strp2.' ...
    strp3 strp4]);
csvwrite(strcat('bat_simulation',strp5,'.csv'),log_b);
csvwrite(strcat('ctl_simulation',strp5,'.csv'),[time_v.' log_c1 log_c2]);


%% functions


function [c2_chain] = mpc(min_c1,max_c1,min_c2,max_c2,c1,c2,N,eu,b0,b,...
                          qc,int_v,q0,Ad,B,C,u,est_u,k,delta_T)

    import casadi.* % import casadi for optimal control

    state = MX.sym('state',7,1); % define the states
                                 % model order is 3 (7 states)
    input = MX.sym('input',2,1); % one input, the computation
                                 % param (we replan the path
                                 % param out of mpc checking
                                 % the battery time, as in MPC 
                                 % check for the horizon, not 
                                 % total time needed to drain 
                                 % battery)
    qq = [state(1,1);state(2,1);state(3,1);state(4,1);...
          state(5,1);state(6,1);state(7,1)]; % model variables
    uu = [input(1,1);input(2,1)]; % control
                
    dq = Ad*qq+B*uu;
                
    fode = Function('fode',{state,input},{dq},...
        {'state','input'},{'d_state'});
    STATE = MX.sym('STATE',size(state)); % for the casadi fun
    INPUT = MX.sym('INPUT',size(input));
    QF = MX.sym('QF',size(state)); % output of the integration
    dt = MX.sym('dt');
                
    k1 = fode(STATE,INPUT);
    k2 = fode(STATE+dt/2*k1,INPUT);
    k3 = fode(STATE+dt/2*k2,INPUT);
    k4 = fode(STATE+dt*k3,INPUT);
    QF = STATE+dt/6*(k1+2*k2+2*k3+k4);
    F_RK4 = Function('F_RK4',{STATE,INPUT,dt},{QF},...
        {'STATE','INPUT','dt'},{'QF'}); % rk4 fix step
                
    interval = k:k+N; % no k+1 in the formula as we add the 
                      % initial condition implicitly
    hh = length(interval);
    dt = delta_T;

                
    opti = casadi.Opti(); % define opt problem
    Q = opti.variable(7,hh);
    U = opti.variable(2,hh-1);
                
    initq = q0;                
    initcon = Q(:,1) == initq;
    battcon = C*Q(:,1) <= b0*qc*int_v;
    opti.subject_to(initcon);
    opti.subject_to(battcon);
                
    opti.subject_to(U(2,1) == c2);
    opti.subject_to(U(1,:) == c1);
    opti.subject_to(min_c2 <= U(2,:) <= max_c2);
    opti.subject_to(min_c1 <= U(1,:) <= max_c1);
                
    con = {initcon};
    btc = {battcon};
                
    eeu = eu; % control estimate (from model)
    eeeu =  est_u(c1,c2);
                
    for jj=1:hh-1
        
        result_rk4 = F_RK4(Q(:,jj),u(eeeu(2),eeu(2)),dt); % model dynamics
        b0 = b0+dt*b(C*result_rk4); % battery dynamics

        if mod(jj-1,1/dt) == 0 % every sum in the MPC
                     
            if jj < hh
                eeeu = est_u(U(1,jj),U(2,jj));
            end
                        
                con{jj+1} = Q(:,jj+1) == F_RK4(result_rk4,...
                    u(eeeu(2),eeu(2)),dt); % state const
                        
                btc{jj+1} = C*Q(:,jj+1) <= b0*qc*int_v; % battery const
                opti.subject_to(con{jj+1});
                opti.subject_to(btc{jj+1});
                        
        end                    
                                       
        eeu = eeeu;
    end
                
    L = sum(U(2,:).^2);
                
    opti.minimize(-L);
    opti.solver('ipopt');
                
    try
        sol = opti.solve();
        c2_chain = sol.value(U); % optimal u on N
        c2_chain = c2_chain(2,:);
    catch
        c2_chain = ones(1,N)*min_c2; % there is no control 
                                     % which sattisfies consts
    end            

end

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
    ddot_pd = (E-ke*e*eye(2))*H*dot_p-ke*n'*dot_p*n; % (10)
    ddot_pdhat = -E*(dot_pd*dot_pd')*E*ddot_pd/norm(dot_pd)^3; % (9)

    dot_Xid = ddot_pdhat'*E*dot_pd/norm(dot_pd); % (13)

    u_theta = dot_Xid+kd*dot_p'*E*dot_pd/(norm(dot_p)*norm(dot_pd)); % (16)

end

function [A,B,C] = build_model(omega,r)
%BUILD_MODEL Creates the simplified model (local implementation from
%build_model)
    m = 2*r+1;

    Aj = @(omega, j) [0 omega*j;-omega*j 0];
    A  = zeros(m);
    A(1,1) = 0;
    B = [0 1]; % two controls
    C = [1];

    for i = 1:r
        A(2*i:2*i+1,2*i:2*i+1) = Aj(omega,i); 
        C = [C 1 0];
        B = [B;0 0;0 0];
    end

end

