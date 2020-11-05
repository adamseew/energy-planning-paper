function [k pos pow est ua] = simulate(varphi, trigger, trajparams, mpcparams, si, ci, animate)
%SIMULATE Simulation wrapper function
%   
% Inputs:
%   varphi:   the set of TEE [strings, format: rotation;ke;equation]
%   trigger:  the set of triggering point that tells when to go from TEE i 
%             to TEE i+1
%   trajparams: parameters about the trajectory:
%     |-vehspeed: vehicle speed [m/s]
%     |-windspeed:wind speed [m/s]
%     |-vehdir:   vehicles direction [angle]
%     |-winddir:  wind direction
%     |-strat:    initial x coordinate [m]
%     |-end:      initial y coordinate [m]
%     |-maxpw:    max power that the vehicle can reach [W]
%     |-minpw:    min power [W]
%     |-trigeps:  the radius of the circle surrounding the triggering point
%     |___________which allows to change the TEE from i to i+1
%   mpcparams: parameters about the model and the algorithm
%     |-order:    order of the Fourier series the model is derived from
%     |-xi:       characteristic time
%     |-epsilon:  treshold which says when to use model, when data for
%     |           the state observer
%     |-horizon:  optimization horizon (how far in the future the
%     |___________algorithm looks
%   si:       the QoS set (must containt all the possible values per each 
%             t)
%   ci:       the TEEs set (same applies as above)
%   animate:  show animation
%
% Outputs:
%   k:        time [sec]
%   pos:      vector of positions [2m]
%   pow:      vector of powers [W]
%   est:      estimated energy with state observer [W]
%   ua:       optimal control input sequence derived with MPC
%

    vehspeed = trajparams(1); 
    windspeed = trajparams(2); 
    vehdir = trajparams(3); 
    winddir = trajparams(4);
    startx = trajparams(5);
    starty = trajparams(6);
    maxpw = trajparams(7);
    minpw = trajparams(8); 
    trigeps = trajparams(9);
    
    r = mpcparams(1);
    xi = mpcparams(2);
    eps = mpcparams(3);
    N = mpcparams(4);
       
    [A, B, C, u, h] = build_model(r, xi, si, ci); 

    % initial guess
    [q0, P0, Q, R] = guess_kf(minpw, size(A, 1));

    % MPC tuning matrices
    Rm = eye(size(B,2));
    Qm = eye(size(q0,1));
    Pf = Qm;

    addpath(genpath('position'));
    addpath(genpath('energy'));

    E = [0 -1; 1 0];; % default rotation matrix
    ke = 0.0003; % default gain to adjust speed of convergence
    
    k = 0;
    i = 1;
    
    triggereps = trigeps; % tollerance to the triggering position
    
    pos = [];
    pow = [];
    
    nowpos = [startx starty];
    pdanglelist = [];
    
    % contribution of the wind. This doesn't change as winspeed and
    % direction is constant
    windx = .1 * windspeed * cosd(winddir);
    windy = .1 * windspeed * sind(winddir);
    
    y = []; % model output
    q = []; % model state

    for varphii = transpose(varphi) % per each TEE
        
        varphii = split(varphii, ";");
        
        fprintf('varphi%d((%d,%d)):=%s\n', i, nowpos(1), nowpos(2), varphii(3));
        
        E = [0 -1; 1 0];
        if contains(varphii(2), '270')
            E = [0 1; -1 0];
            disp('rotation=270');
        else
            disp('rotation=90');
        end
        
        ke = str2double(varphii(1));
        
        fprintf('ke=%f\n', ke);
        
        while true

            pos = [pos; nowpos];
            
            % vector field
            [dpd, pdangle] = build_gdn2(E, ke, str2sym(varphii(3)), nowpos); 
            
            % dpd is the ideal direction... We have wind though
            % now if the time is sampled every second, the below expression
            % actually indicates the offset from the original location
            
            posx = .1 * vehspeed * cosd(pdangle);
            posy = .1 * vehspeed * sind(pdangle);
            
            % new position
            nowpos = [nowpos(1) + windx + posx, nowpos(2) + windy + posy]; 
            
            % reached the triggering point, going to TEE i+1
            if all(abs(nowpos - trigger(i,:)) <= triggereps) 
                break;
            end
            
            % plotting the path
            subplot(2,2,[1 2]);
            plot(pos(:, 1), pos(:, 2), 'Color', 'r', 'LineWidth', 1.2)
            
            
            % produce a simulated energy value
            % first, let's calculate the angle between the two vectors
            angle = acosd(...
                (posx * windx + posy * windy) / ...
                (sqrt(posx^2 + posy^2) * sqrt(windx^2 + windy^2))...
                         );
            % angle is nan? no windspeed then; so the contribution of the
            % wind is null
            if isnan(angle)
                angle = 0;
            end
            
            % if the angle is 0, lowest energy (tail wind). 180, heighest
            % (head wind)
            simpow = minpw + (maxpw - minpw) * .5 * abs(1 - cosd(angle)) * ... 
                windspeed / vehspeed; % winddir and windspeed effect
            
            % before it was
            %simpow = maxpw - (.5 * (maxpw - minpw) + ...
            %    .5 * (maxpw - minpw) * cosd(winddir - pdangle)); % wind
            
            pow = [pow; simpow];
            
            % plotting the energy and the estimated energy
            subplot(2,2,3);
            plot(linspace(0, k/10, size(pow, 1)), pow(:, 1))
            
            % estimating the energy every second
            if and(k ~= 0, mod(k, 10) == 0)

            
                % output MPC
                %[~, ~, ~, yy, P1, q0] = output_mpc(A, B, C, si, q0, P0, ...
                %    Q, R, pow, ...
                %    Rm, Qm, Pf, k/10, eps, N);
                
                %y = []; % instantaneous energy consumption (evolution in 
                        % time)
                %q = []; % state (also evolution in time)
                %ua = []; % control input sequence
                
                % iteration over all the possible controls
                % this instruction corresponds to argmax in algorithm (line
                % 3)
                %v = [];
                %for current_si = transpose(si(any(si(:, 1) == i, 2), end))
                %    qk = q0;
                    
                %    vv = .5 * (trasnpose(qk) * Qm * qk + transpose(current_si) * Rm * current_si);
                    
                    % this instruction corresponds to the sum 
                %    for l = 1 : N - 1
                        
                %    end
                %end
                    
                % first we need to interate over all the ppp
            
                % plotting the estimated energy
                %subplot(2,2,4);
                %plot(linspace(0, k/10, size(yy, 1)), yy(:, 1))
                
                % for now I do just the estimation with the model and later
                % implement the rest of the MPC
                yy = [];
                P1 = [];
                
                [yy, qq0, P1] = estimate_kf(A, B, C, u, q0, P0, Q, R, ...
                    pow(k), k, eps);
                
                y = [y; yy];
                q0 = qq0;
        
                if ~isempty(P1) % if it's empty, no KF just evolutin (<= epsilon)
                    P0 = P1;
                    fprintf("kf at time %d\n", k/10);
                else
                    fprintf("model evolution at time %d\n", k/10);
                end
                q = [q; q0];
                
                subplot(2,2,4);
                plot(linspace(0, k/10, size(y, 1)), y(:, 1))
            end
            
            pause(.1)
            
            pdanglelist = [pdanglelist; pdangle];
            
            k = k + 1;
        end
        
        i = i + 1;
    end

end

