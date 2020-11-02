%OP5
% Simulation of the Energy-Aware Dynamic Mission Planning Algorithm
%


%% Build the model 

fprintf(['[ OP5 ] \n' ...
         '        \n']);
disp( '[     ] Build the model');
fprintf(['[   ! ] Dependencies: mission specification, the value of the\n' ...
         '        model from the modeling tool, t\n']);

% The former computational model (per each time step all the
% possibilities)
ckgck = build_ck(file3, file2, t);

clear file2 file3;

gck = [];
for i = 1:t(end)
    cckgck = ckgck(any(ckgck(:, 1) == t(i), 2), end);
    
    gck = [gck; cckgck(end)];
end

%TODO: the former mechanical components of the control
mk = zeros(size(t, 2), 1);

%% Output MPC

disp('[     ] Output MPC');
disp('[   ! ] Dependencies: A, B, C, control vector u, sensor data meas, t');

eps = input('[   ? ] Input epsilon: ');
if isempty(eps)
    eps = 1/3; % default eps [W tolerance]
end

N = input('[   ? ] Input the horizon: ');
if isempty(N)
    N = 10; % default horizon
end

% gck is the control input sequnce (OP5 adapts the computations)
    
[A, B, C, u, h] = build_model(r, xi, gck, mk); 

% initial guess
[q0, P0, Q, R] = guess_kf(meas(1), size(A, 1));

% MPC tuning matrices
Rm = eye(size(B,2));
Qm = eye(size(q0,1));
Pf = Qm;

% output MPC
[y, q, ua, meas] = output_mpc(A, B, C, ckgck, q0, P0, Q, R, meas, ...
    Rm, Qm, Pf, t, eps, N);
