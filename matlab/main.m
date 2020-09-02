
%main.m
% Simulation in Matlab, use data directory to pass the required data in
% the simulation. If you use the combo c_e (computational_energy.cvs...), 
% m_e, m_s in simulation1, then use r = 3, column = 1, and ts = 200
%

%%
disp('Model simulation')

r = input('Fourier series order r: ');

% getting mechanical energy data

disp('Input mechanical energy data [csv]');
[bn, folder] = uigetfile('.csv');
file = csvread(fullfile([folder, bn]));

column = input('Which column contains energy: ');
file = file(:, [column]);

ts = input('At what time step [ms]: ');
tf = (1 / ts) * 1000 * size(file, 1);
file = imresize(file, [100, 1]);
t = transpose(linspace (0, tf, 100));

clear column;

% getting computational energy data
disp('Input computational energy data (from powprof) [csv]');
[bn, folder] = uigetfile('.csv');
file2 = csvread(fullfile([folder, bn]));

disp('Input mission specification [csv]');
[bn, folder] = uigetfile('.csv');
file3 = csvread(fullfile([folder, bn]));

% computational model
ckgck = build_ck(file3, file2, t);

% initial guess
j = 2*r + 1;
q0 = ones(j, 1);
P0 = eye(j);

% process noise and sensor noise
Q = eye(j);
R = 1;

% TODO: this selects just the highest possible control for each time
% interval, which is not optimal; this would be part of the optimal control
% policy once it's done
gck = [];
for i = 1:size(t, 1)
    cckgck = ckgck(any(ckgck(:, 1) == t(i), 2), end);
    
    gck = [gck; cckgck(end)];
end

%TODO: the mechanical components of the control
mk = zeros(size(t, 1), 1);

% predicted energy evolution
[y, qq0] = build_model(r, file,  gck, mk, q0, P0, Q, R, t, 0.001);

%%
plot(t, file + gck);
clear file file2 file3;
hold on;
plot (t, y);

legend('data', 'observer');

%%

syms x y;
varphi = ((x - 200) * cosd(12) + (y - 200) * sind(12))^2 / 120^2 + ((x - 200)*sind(12) - (y - 200) * cosd(12))^2 / 170^2 - 1;

min = 0;
max = 400;

% create all the points I want to use to build gdn
[xp yp] = meshgrid(min:(abs(min) + abs(max)) / 50:max, min:(abs(min) + abs(max)) / 50:max);
points = reshape(cat(2, xp', yp'), [], 2);
clear xp yp;

E = [0 1; -1 0];
ke = 1.5;
[dpd, pdangle] = build_gdn2(E, ke, varphi, points);
plot_gdn2(E, ke, varphi, points, pdangle, 0, 400, 10000, 20, [400; 400])




