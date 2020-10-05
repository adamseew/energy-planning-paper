function [A, B, C, u, h] = build_model(r, xi, gck, mk)
%BUILD_MODEL.M Builds the energy model using Kalman Filter as state 
%   observer
%
% Inputs:
%   r:  the Fourier series order [Equation (6)]
%   xi: characteristics time
%   gck:control the value of the power interrogating powprof with the 
%       QoS parameters
%   mk: control TEEs parameters values
%
% Outputs:
%   A:  state matrix
%   B:  control matrix
%   C:  output matrix
%   u:  control vector
%   h:  Fourier series
%

    % system definition
    j = 2*r + 1;

    An = @(n) [ 0 n/xi ; -n^2/xi^2 0 ];
    A  = zeros(j);
    A(1,1) = 1;
    C = [1];
    qm = zeros(1, j);
    syms q [j 1];
    syms k;
    
    h = q(1); 
    
    for i = 1:r
        A(2 * i : 2 * i + 1, 2 * i : 2 * i + 1) = An(i); 
        C = [C 1 0];
        h = h + q(2 * i) * cos(i * k / xi) + ...
                q(2 * i + 1) * sin(i * k / xi);
    end
    
    clear k;
    
    % control action
    u = @(k) [gck(k); transpose(mk(k))];
    B = [ 1                 zeros(1, size(mk, 2)); ...
          zeros(j - 1, 1)   zeros(j - 1, size(mk, 2))  ];

end
