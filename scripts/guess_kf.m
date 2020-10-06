function [q0, P0, Q, R] = guess_kf(meas1, j)
%GUESS_KF Initial guess of the parameters for KF (and AKF)
%
% Inputs:
%   meas1: the first measurement from the sensor
%   j:     desired size of the state
%
% Outputs:
%   q0:    initial state guess
%   P0:    variance of the guessed estimate
%   Q:     uncertainty error covariance guess
%   R:     measurement error covariance guess
%

    q0 = ones(j, 1);

    % following operation put a lot of weight on the first element, such as
    % in the original series
    
    q0(1) = meas1 * ((j -1) / j);
    q0(2:end) = meas1 * (1 / ((j - 1) * j));
    
    P0 = eye(j) .* q0 * 10^-1;
    
    % process noise and sensor noise
    Q = eye(j) .* q0 * 10^-1;
    R = 10 * meas1;    

end

