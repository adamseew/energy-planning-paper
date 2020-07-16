
%% system definition

r = input('Fourier series order r: ');
j = 2*r + 1;

An = @(n) [ 0 1 ; -1 0 ];
A  = zeros(j);
A(1,1) = 1;
C = [1];

for i = 1:r
    A(2 * i : 2 * i + 1, 2 * i : 2 * i + 1) = An(i); 
    C = [C 1 0];
end

%% getting energy data

file = csvread(uigetfile('.csv'));
column = input('Which column contains energy: ');
file = file(:, [column]);
ts = input('At what time step [ms]: ');
tf = (1 / ts) * 1000 * size(file, 1);
file = imresize(file, [100, 1]);
t = transpose(linspace (0, tf, 100));

%% initial guess

q0 = ones(j, 1);
P0 = eye(j);
%u=...

% process noise and sensor noise
Q = eye(j);
R = 1;

% predicted energy evolution
y = [];

%% kalman filter estimatation
 
for y0_sensor = transpose(file)
    
    % prediction
    q1_minus = A * q0; %+B*u0
    P1_minus = A * P0 * transpose(A) + Q;
    
    % estimation
    K1 = (P1_minus * transpose(C)) / (C * P1_minus * transpose(C) + R);
    q1 = q1_minus + K1 * (y0_sensor - C * q1_minus);
    P1 = (eye(j) + K1 * C) * P1_minus;
    y0_estimate = C * q1;
    
    % updating values for the next iteration
    q0 = q1;
    P0 = P1;
    y = [y y0_estimate];
end

%% plots

plot(t, file);
hold on;
plot (t, y);

legend('data', 'observer');



