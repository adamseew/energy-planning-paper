
function [u, g] = read_spec(spec, prof, t) 
%mission.m Reads a mission specification file along with the output from
%   powprof. Thus builds a vector of all the possible controls and their 
%   power consumption. The control vector is normalized over the systems
%   time, for the purpose of simulation.
%
% Inputs:
%   spec   : mission specification
%   prof   : output from powprof (I expect this comes unchanged from
%            powprof [which means ordered])
%   t      : time used for the purpose of simulation
%
% Outputs:
%   comp   : The computational energy component u = {ui | ui = QoS for 
%                                                        all i in {0, ..., 
%                                                                  k}     }
   
    u = [];
    
    for i = 1:size(spec, 1)
        % iterates for each line in the mission specification
        
        for j = spec(i, 1):spec(i, 2)
            % iterates for each time interval
            
            u1 = [];
            for k = 3:2:size(spec, 2)
                % makes all the possible combinations of controls of QoS
                
                uu = spec(i, k): ...
                     prof(2, end - 1) - prof(1, end - 1):spec(i, k + 1);
                % for the purpose of simulation the same step as in powprof
                % is used (which can be easily adapted when using powprof)
                
                if (isempty(u1))
                    u1 = transpose(uu);
                else
                    u1 = repmat(u1, size(uu, 2), 1);
                    uu = repmat(uu, size(u1, 1) / size(uu, 2), 1);
                    uu = uu(:);
                    u1 = [ u1 uu ];
                end
            end
            
            tt = ones(size(u1, 1), 1) * j;
            
            u1 = [ u1 zeros(size(u1, 1), 1) ];
            % filling the control with the energy data
            for k = 1:size(prof, 1)                
                u1(all( ...
                        u1(:, 1:end - 1) == prof(k, 1:end - 1), 2 ...
                      ), end) = prof(k, end);
            end
            
            u = [ u; tt u1 ];
            
            clear tt uu u1 ui;
        end
    end
    
end