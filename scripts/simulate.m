function [k pos pow] = simulate(varphi, trigger, vehspeed, windspeed, ...
    vehdir, winddir, startx, starty, maxpw, minpw, animate)
%SIMULATE Simulation wrapper function
%   
% Inputs:
%   varphi:   the set of TEE [strings]
%   trigger:  the set of triggering point that tells when to go from TEE i to
%             TEE i+1
%   vehspeed: vehicle speed [m/s]
%   windspeed:wind speed [m/s]
%   vehdir:   vehicles direction [angle]
%   winddir:  wind direction
%   strat:    initial x coordinate [m]
%   end:      initial y coordinate [m]
%   maxpw:    max power that the vehicle can reach [W]
%   minpw:    min power [W]
%   animate:  show animation
%
% Outputs:
%   k:        time [sec]
%   pos:      vector of positions [2m]
%   pow:      vector of powers [W]
%

    addpath(genpath('position'));
    addpath(genpath('energy'));

    E = [0 -1; 1 0];; % default rotation matrix
    ke = 0.0005; % default gain to adjust speed of convergence
    
    k = 0;
    i = 1;
    
    triggereps = 10; % tollerance to the triggering position
    
    pos = [];
    pow = [];
    
    olddpd = 0;
    
    nowpos = [startx starty];
        
    for varphii = transpose(varphi) % per each TEE
        
        while true

            pos = [pos; nowpos];
            
            % vector field
            [dpd, pdangle] = build_gdn2(E, ke, str2sym(varphii), nowpos); 
            
            % dpd is the ideal direction... We have wind though
            % now if the time is sampled every second, the below expression
            % actually indicates the offset from the original location
            windx = windspeed * cosd(winddir);
            windy = windspeed * sind(winddir);
            posx = vehspeed * cosd(pdangle);
            posy = vehspeed * sind(pdangle);    
    
            % new position
            nowpos = [nowpos(1) + windx + posx, nowpos(2) + windy + posy]; 
            
            % reached the triggering point, going to TEE i+1
            if all(abs(nowpos - trigger(i,:)) <= triggereps) 
                break;
            end
            
            % produce a simulated energy value; if the drone doesn't turn,
            % no energy is needed due to turn (minus first contribution)
            % if there is tailwind, no energy needed due to wind (minus
            % second)
            simpow = maxpw - ...
                .5 * (maxpw - minpw) * abs(pdangle - olddpd) / 180 - ...
                .5 * (maxpw - minpw) * abs(windspeed - pdangle) / 180;
            
            pow = [pow; simpow];
            
            k = k + 1;
        end
        
        i = i + 1;
    end

end

