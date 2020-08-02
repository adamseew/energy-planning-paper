
function [ckgck] = build_ck(spec, prof, t) 
%build_ck.m Reads a mission specification file along with the output from
%   powprof. Builds a vector of all the possible computational controls 
%   [Ck in Equation (12)] and their  power consumption [g(Ck) in 
%   Equation (13)]. The control vector is normalized over the systems
%   time, for the purpose of simulation.
%
% Inputs:
%   spec   : mission specification
%   prof   : output from powprof (I expect this comes unchanged from
%            powprof [which means ordered])
%   t      : time used for the purpose of simulation
%
% Outputs:
%   ckgck  : All the possible combinations of the QoS Ck, along with their
%            instantaneous power consumption from powprof g(Ck)
%   
    u = [];
    
    % iterates for each line in the mission specification
    for i = 1:size(spec, 1)
        
        % iterates for each time interval
        for j = spec(i, 1):spec(i, 2)
            
            % just for the purpose of simulation; you generate all the
            % possible combinations of computational control for the
            % simulated time (and not all the times) 
            if (~any(floor(t) == j))
                continue;
            end
            jj = t(any(floor(t) == j, 2));
            jj = jj(1);
            
            % makes all the possible combinations of controls of QoS
            u1 = [];
            for k = 3:2:size(spec, 2)
                
                % for the purpose of simulation the same step as in powprof
                % is used (which can be easily adapted when using powprof)
                uu = spec(i, k): ...
                     prof(2, end - 1) - prof(1, end - 1):spec(i, k + 1);
                
                if (isempty(u1))
                    u1 = transpose(uu);
                else
                    u1 = repmat(u1, size(uu, 2), 1);
                    uu = repmat(uu, size(u1, 1) / size(uu, 2), 1);
                    uu = uu(:);
                    u1 = [ u1 uu ];
                end
            end
            
            % filling the control with the energy data
            u1 = [ u1 zeros(size(u1, 1), 1) ];
            
            for k = 1:size(prof, 1)                
                u1(all( ...
                        u1(:, 1:end - 1) == prof(k, 1:end - 1), 2 ...
                      ), end) = prof(k, end);
            end
            
            tt = ones(size(u1, 1), 1) * jj;
            
            u = [ u; tt u1 ];
            
            clear tt uu u1 ui;
        end
    end
    
    ckgck = u;
    
end
