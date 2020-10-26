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

% why meas + gck? Including both contributions; the former conputational
% and mechanical energy from the sensor

cckgck = [];

% MPC tuning matrices
Rm = eye(size(B,2));
Qm = eye(size(q0,1));
Pf = Qm;

y = [];
q = [];
ua = [];
for k=1:t(end)
   
    % line 28
    cckgck = ckgck(any(ckgck(:, 1) == k, 2), end);
        
    % iteration over all the possible controls
    v = [];
    for current_gck=transpose(cckgck)
        vv = 0;    
        qk = q0;
        
        vv = vv + .5*(transpose(qk)*Qm*qk+current_gck*current_gck);
        
        for l=1:N-1
            for all_gck=transpose(ckgck(any(ckgck(:, 1) == k + l, 2), end))
                all_gck = [all_gck; 0]; % simulate also the adaptations
                vv = vv + .5*(transpose(qk)*Qm*qk+transpose(all_gck)*Rm*all_gck);
            end
        
            [yk qk] = evolve_sys(A, B, C, current_gck, qk, 1);
            qk = qk(:,end); % keep just last one, others not necessary
        end
        
        vv = vv + .5*(transpose(qk)*Pf*qk);
        
        % this is a set of ordered couples, second one is the value of the
        % cost function (see line 3), the first one is the control
        v = [v; vv current_gck]; % just the computations for now
    end
    
    [vv i] = max(v(:, 1));
    maxu = v(i, 2);
    meas(k) = meas(k) + maxu; % adding the former computational energy
    ua = [ua maxu];
    if size(ua,2) >= 2
        maxu = @(k) [maxu - ua(1,end - 1); 0];
    else
        maxu = @(k) [maxu; 0]; % a silly workaround due to the old estimator
    end
    yy = [];
    P1 = [];
    [yy, qq0, P1] = estimate_kf(A, B, C, maxu, q0, P0, Q, R, meas(k), k, eps);
    y = [y; yy];
    q0 = qq0;
    if ~isempty(P1)
        P0 = P1;
    end
    q = [q; q0];
end