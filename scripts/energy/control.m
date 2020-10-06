function [u] = control(k, gck, mk)
%CONTROL Defines the control action (wrapper)
%
% Inputs:
%   gck:control the value of the power interrogating powprof with the 
%       QoS parameters
%   mk: control TEEs parameters values
%
% Outputs:
%   u:  control vector
%

    if (k ~= 1)
        u = [gck(k) - gck(k - 1); transpose(mk(k)) - transpose(mk(k - 1))];
    else
        u = [gck(k); transpose(mk(k))];
    end
end

