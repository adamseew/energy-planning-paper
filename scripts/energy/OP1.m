%OP1
% Simulation of the Energy-Aware Dynamic Mission Planning Algorithm
%
% Non animated simulation / fixed case: with no TEEs controls, and the 
% highest possible QoS controls at each time step


%% Build the model 

fprintf(['[ OP1 ] Non animated simulation / fixed case: with no TEEs\n' ...
         '        controls, and the highest possible QoS controls at each\n' ... 
         '        time step\n']);
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

disp('[     ] Estimate the state');
disp('[   ! ] Dependencies: A, B, C, control vector u, sensor data meas, t');

% initial guess
[q0, P0, Q, R] = guess_kf(meas(1), size(A, 1));

% why meas + gck? Including both contributions; the former conputational
% and mechanical energy from the sensor
meas = meas + gck;

[y, q] = estimate_kf(A, B, C, u, q0, P0, Q, R, meas, t, 0);

