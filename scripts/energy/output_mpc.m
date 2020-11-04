function [y, q, ua, meas, P1, q0] = output_mpc(A, B, C, ck, q0, P0, Q, R, meas, ...
    Rm, Qm, Pf, t, eps, N)
%OUTPUT_MPC Output MPC algorithm (see MPC Rawlings and Mayne, 2009, p. 382)
%
% Inputs:
%   A:   state matrix
%   B:   control matrix
%   C:   output matrix
%   ck:  the QoS set (must containt all the possible values per each t)
%   q0:  initial guess of the estimate of the state
%   P0:  variance of the guessed estimate
%   Q:   uncertainty error covariance
%   R:   measurement error covariance
%   meas:data from the sensor (i.e., measured instantaneous energy
%        consumption, the former mechanical energy from sensor)
%   Rm:  MPC tuning param (output)
%   Qm:  MPC tuning param (state)
%   Pf:  MPC tuning param (final step ouput)
%   t:   time vector
%   eps: epsilon (if 0 then normal KF is used instead of adaptative KF)
%   N:   MPC horizon
%
% Outputs:
%   y:   instantaneous energy consumption (evolution in time)
%   q:   state (also evolution in time)
%   ua:  control input sequence (optimal)
%   meas:updated data from sensor with the computational energy
%   P1:  next variance of the guessed estimate
%   q0:  next state
%

    y = []; % output
    q = []; % state

    ua = []; % control input

    for k = 1:t(end)
        
        % iteration over all the possible controls
        v = [];
        for current_gck = transpose(ck(any(ck(:, 1) == k, 2), end))
            vv = 0;    
            qk = q0;
        
            vv = vv + .5*(transpose(qk)*Qm*qk+current_gck*current_gck);
        
            for l=1:N-1
                for all_gck = transpose(ck(any(ck(:, 1) == k + l, 2), end))
                    all_gck = [all_gck; 0]; % simulate also the adaptations
                    vv = vv + .5 * (transpose(qk) * Qm * qk + ...
                        transpose(all_gck) * Rm * all_gck); % step function
                end
        
                [yk qk] = evolve_sys(A, B, C, current_gck, qk, 1);
                qk = qk(:, end); % keep just last one, others not necessary
            end
        
            vv = vv + .5 * (transpose(qk) * Pf * qk); % final step function
        
            % this is a set of ordered couples, second one is the value of 
            % the cost function (see line 3), the first one is the control
            v = [v; vv current_gck]; % just the computations for now
        end
    
        [vv i] = max(v(:, 1)); % maxof the set (max arg). Ideally you'd
                               % implement KKT to solve QP
        maxu = v(i, 2);
        meas(k) = meas(k) + maxu; % adding the former computational energy
        ua = [ua maxu];
        
        if size(ua, 2) >= 2
            maxu = @(k) [maxu - ua(1,end - 1); 0];
        else
            maxu = @(k) [maxu; 0]; % a workaround due to the old estimator
        end
        
        yy = [];
        P1 = [];
        
        [yy, qq0, P1] = estimate_kf(A, B, C, maxu, q0, P0, Q, R, ...
            meas(k), k, eps);
        y = [y; yy];
        q0 = qq0;
        
        if ~isempty(P1) % if it's empty, no KF just evolutin (<= epsilon)
            P0 = P1;
        end
        q = [q; q0];
    end

end


