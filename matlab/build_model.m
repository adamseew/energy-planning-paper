
function [y, qq0] = build_model(r, mech, gck, mk, q0, P0, Q, R, t, eps)
%build_model.m Builds the energy model [Equation (9b)] using Kalman Filter 
%   as state observer [Equations (10--11)]
%
% Inputs:
%   r      : the Fourier series order [Eqaution (6)]
%   mech   : mechanical energy vector [W]. I.e., the data from flying the
%            drone in simulation or data from real flights logging the
%            flight controller
%   u      : control [Equation (9b)]
%   q0     : initial guess of the estimate of the state
%   P0     : vatiance of the guessed estimate
%   Q      : process noise covariance
%   R      : sensor noise covariance
%   t      : time used for the purpose of simulation
%   eps    : used to return qq0
%
% Outputs:
%   y      : energy evolution
%   qq0    : the first estimated state value [Equation (8)] differs from
%            mech < eps
%

    % TODO: u!
    % system defintion
    j = 2*r + 1;
    xi = 10;

    An = @(n) [ 0 n/xi ; -n^2/xi^2 0 ];
    A  = zeros(j);
    A(1,1) = 1;
    C = [1];

    for i = 1:r
        A(2 * i : 2 * i + 1, 2 * i : 2 * i + 1) = An(i); 
        C = [C 1 0];
    end
    
    % control action [Equation (17)]
    u = @(k) [gck(k); transpose(mk(k))];
    B = [ 1                 zeros(1, size(mk, 2)); ...
          zeros(j - 1, 1)   zeros(j - 1, size(mk, 2))  ];

    clear i, An;
    
    % predicted energy evolution
    y = [];

    qq0 = [];
    
    % kalman filter estimatation
    k = 1;
    for y0_sensor = transpose(mech)
        % getting the control at time k
        u0 = u(k);
        
        % prediction
        q1_minus = A * q0 + B*u0;
        P1_minus = A * P0 * transpose(A) + Q;
    
        % estimation
        K1 = (P1_minus * transpose(C)) / (C * P1_minus * transpose(C) + R);
        q1 = q1_minus + K1 * (y0_sensor + u0(1) - C * q1_minus);
        P1 = (eye(j) + K1 * C) * P1_minus;
        y0_estimate = C * q1;
    
        % updating values for the next iteration
        q0 = q1;
        P0 = P1;
        k = k + 1;
        
        y = [y y0_estimate];
        
        if (and(isempty(qq0), abs(y0_estimate - y0_sensor) < 0.001))
            qq0 = q0;
        end
    end

    clear j showed;

end
