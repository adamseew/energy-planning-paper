function [y, q] = evolve_sys(A, B, C, u, q0, t)
%ESTIMATE_KF Discrete-time linear Kalman filter state estimation
% (see: Dan Simon, Optimal State Estimation [ISBN: 9780471708582], p.128)
%
% Inputs:
%   A:   state matrix
%   B:   control matrix
%   C:   output matrix
%   u:   control vector (must contain one value per each item in t)
%   t:   time vector
%
% Outputs:
%   y:   instantaneous energy consumption (evolution in time)
%   q:   state (also evolution in time)
%

    % predicted instantaneous energy consumption evolution in time
    y = [];
    
    % estimated state evolution in time
    q = [];
    
    for k = t
        % getting the control at time k
        u0 = u(k);
        
        % prediction
        q1 = A * q0 + B * u0;
        y1 = C * q1;
        
        y = [y y1];
        q = [q q1];
        
        q0 = q1;
        
    end

end

