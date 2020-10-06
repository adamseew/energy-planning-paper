function [y, q] = estimate_kf(A, B, C, u, q0, P0, Q, R, meas, t, eps)
%ESTIMATE_KF.M Discrete-time linear Kalman filter state estimation
% (see: Dan Simon, Optimal State Estimation [ISBN: 9780471708582], p.128)
%
% Inputs:
%   A:   state matrix
%   B:   control matrix
%   C:   output matrix
%   u:   control vector (must contain one value per each item in t)
%   q0:  initial guess of the estimate of the state
%   P0:  variance of the guessed estimate
%   Q:   uncertainty error covariance
%   R:   measurement error covariance
%   meas:data from the sensor (i.e., measured instantaneous energy
%        consumption, the former mechanical energy from sensor)
%   t:   time vector
%   eps: epsilon (if 0 then normal KF is used instead of adaptative KF)
%
% Outputs:
%   y:   instantaneous energy consumption (evolution in time)
%   q:   state (also evolution in time)
%
    
    % predicted instantaneous energy consumption evolution in time
    y = [];
    
    % estimated state evolution in time
    q = [];
    
    ddd = [];
    
    % Discrete-time KF estimation
    k = 1;
    for y0_sensor = transpose(meas)
        % getting the control at time k
        u0 = u(k);
        
        % prediction
        q1_minus = A * q0; %+ B * u0;
        
        if and(eps > 0, ...
               abs(y0_sensor - C * q1_minus) <= eps ... %+ u0(1)) - C * q1_minus) <= eps ...
               ) % use regular system evolution 
           
            q0 = q1_minus;
            y0_estimate = C * q0;
            
            ddd = [ddd; 1 0];
        else % else use KF
                    
            P1_minus = A * P0 * transpose(A) + Q;
    
            % estimation
            K1 = (P1_minus * transpose(C)) / ...
                (C * P1_minus * transpose(C) + R);
            q1 = q1_minus + K1 * ((y0_sensor - C * q1_minus)); %+ u0(1)) - C * q1_minus);
            P1 = (eye(size(q0, 1)) + K1 * C) * P1_minus;
        
            y0_estimate = C * q1;
    
            % updating values for the next iteration
            q0 = q1;
        
            P0 = P1;
            
            ddd = [ddd; 0 1];
        end
        
        k = k + 1;
        
        y = [y y0_estimate];
        q = [q q0];
        
    end

end
