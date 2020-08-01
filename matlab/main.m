
%mainai.m
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

ckgck = build_ck(file3, file2, t);

% initial guess
j = 2*r + 1;
q0 = ones(j, 1);
P0 = eye(j);

% process noise and sensor noise
Q = eye(j);
R = 1;

% predicted energy evolution
[y, qq0] = build_model(r, file,  ckgck(:, end), q0, P0, Q, R, t, 0.001);

%%
plot(t, file);
clear file;
hold on;
plot (t, y);

legend('data', 'observer');
