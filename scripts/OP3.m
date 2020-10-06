%OP3
% Simulation of the Energy-Aware Dynamic Mission Planning Algorithm
%
% Simulation of adaptive KF. KF is activated only when the difference 
% between the output and measurement is greater or equal to epsilon rest 
% same as OP1


%% Build the model 

fprintf(['[ OP3 ] Simulation of adaptive KF. KF is activated only when\n' ...
         '        the difference between the output and measurement\n' ...
         '        is greater or equal to epsilon rest same as OP1\n']);
disp( '[     ] Build the model');
fprintf(['[   ! ] Dependencies: mission specification, the value of the\n' ...
         '        model from the modeling tool, t\n']);

% The former computational model (per each time step all the
% possibilities)
ckgck = build_ck(file3, file2, t);

clear file2 file3;

% TODO: this selects just the highest possible control for each time
% interval, which is not optimal; this would be part of the optimal control
% policy once it's done
gck = [];
for i = t
    cckgck = ckgck(any(ckgck(:, 1) == t(i), 2), end);
    
    gck = [gck; cckgck(end)];
end

clear cckgck ckgck i;

%TODO: the former mechanical components of the control
mk = zeros(size(t, 2), 1);

% Finally builds the model (see eq:state-perf)
% h contains the explicit Fourier series
[A, B, C, u, h] = build_model(r, xi, gck, mk);

clear mk;


%% Estimate the state

disp('[     ] Estimate the state adaptively');
disp('[   ! ] Dependencies: A, B, C, control vector u, sensor data meas, t');

eps = input('[   ? ] Input epsilon: ');
if isempty(eps)
    eps = 1/3; % default eps [W tolerance]
end

% initial guess
[q0, P0, Q, R] = guess_kf(meas(1), size(A, 1));

% why meas + gck? Including both contributions; the former conputational
% and mechanical energy from the sensor
meas = meas + gck;

[y, q] = estimate_kf(A, B, C, u, q0, P0, Q, R, meas, t, eps);

