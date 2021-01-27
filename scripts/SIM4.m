
% simulation from the equations from Hector

clear all; clc; close all;

%% init
delta_T = 1e-3;
final_time = 240;

m = 1; % mass of the aircraft [kg]
W = m*9.8; % weight force [N]

cl = 9.8/15^2; % coefficient to be determined experimentally
cth = 15/50; % same

kp = 5; % positive gain constant to be also determined experimentally
kvv = 5;
th_nominal = 50; % nominal vlalue of the throttle
th_delta = 0;
hd = 25; % desired altitude [m]
u_theta = 0.1;

w = [1;0]; % wind vector (x, y axis velocity

vv = 0; % vertical velocity
sh = 15; % horizontal speed
h = 24; % altitude
p = [0; 0]; % position
theta = 0; % initial angle
vh = sh * [cos(theta);sin(theta)]; % initial velocity
wbx = dot(w,[cos(theta), sin(theta)]);
vs = cth*(th_nominal + th_delta) + wbx; % initial airspeed


log_time = 0:delta_T:final_time;
log_vv = zeros(length(log_time), 1);
log_h = log_vv;
log_theta = log_vv;
log_sh = log_vv;
log_th_delta = log_vv;
log_p = zeros(length(log_time), 2);
log_vh = log_p;

i = 1;
log_vv(i) = vv;
log_h(i) = h;
log_theta(i) = theta;
log_sh(i) = sh;
log_th_delta(i) = th_delta;
log_p(i, :) = p;
log_vh(i, :) = vh;


%% sim

for time = delta_T:delta_T:final_time

% Vertical dynamics
th_delta = kp*(hd - h) - kvv *vv;
wbx = dot(w,[cos(theta), sin(theta)]);
vs = cth*(th_nominal + th_delta) + wbx;
L = cl * vs*vs;
av = 1/m * (L - W);

vv = vv + av*delta_T;
h = h + vv*delta_T;

% Horizontal unicycle model
sh = cth*(th_nominal + th_delta);
pdot = sh*[cos(theta);sin(theta)] + w;

p = p + pdot*delta_T;
theta = theta + u_theta*delta_T;

% Log
i = i+1;
log_vv(i) = vv;
log_h(i) = h;
log_theta(i) = theta;
log_sh(i) = sh;
log_th_delta(i) = th_delta;
log_p(i, :) = p;
log_vh(i, :) = vh;
    
end

%% plots

subplot(2,2,1);
plot(log_p(:, 1), log_p(:, 2))
title('position');
subplot(2,2,2);
plot(log_time, th_nominal + log_th_delta)
title('Throttle');
subplot(2,2,3);
plot(log_time, log_vv);
title('vertical speed');
subplot(2,2,4);
plot(log_time, log_h);
title('altitude');
