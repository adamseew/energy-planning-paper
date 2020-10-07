
function [p, w, segments, circles] = load_pprz(fn, ah)
%LOAD_PPRZ loads the position, power, and the mission (a set of straight 
% lines initial and final points and circles centers and radiuses) from a
% paparazzi autopilot logs
% (see: http://wiki.paparazziuav.org/wiki/Logs)
%
% Inputs:
%   fn:  file name where the log is stored [data]
%   ah:  battery capacity [ah]
%
% Outputs:
%   p:       positions (time [s], relative to home x [m], y [m], 
%            altitude [m])
%   w:       power (time [s], instantaneous power drain [W]
%   segments:straight lines in the mission (time [s], first point x, y [m],
%            second point x, y [m])
%   circles: circles in the mission (time [s], center x, y [m], radius [m])
%

    fid = fopen(fn,'r');
   
    i = 1;
    b_ = 1;
    s_ = 1;
    c_ = 1;
    n_ = 1;

    while feof(fid) == 0
        tline = fgetl(fid);
   
        if  (length(findstr(tline, 'GPS_SOL'))) > 0
            continue;
        end
        
        if ( length(findstr(tline, 'GPS')) ) > 0
            GPS(i, :) = textscan(tline, ...
                '%f %*d %*s %d %d %d %d %d %d %d %d %d %d', 1);
            i = i + 1;
        end
        
        if  (length(findstr(tline, 'NAVIGATION_REF'))) > 0
            continue;
        end
        
        if ( length(findstr(tline, 'NAVIGATION')) ) > 0
            NAV(n_, :) = textscan(tline, ...
                '%f %*d %*s %d %d %d %d %d %d %d %d', 1);
            n_ = n_ + 1;
        end
   
        if ( length(findstr(tline, 'BAT')) ) > 0
            BAT(b_, :) = textscan(tline, ...
                '%f %*d %*s %d %d %d %d %d %d %d %d', 1);            
            b_ = b_ + 1;
        end
        
        if ( length(findstr(tline, 'CIRCLE')) ) > 0
            CIR(c_, :) = textscan(tline, '%f %*d %*s %d %d %d', 1);
            c_ = c_ + 1;
        end
        
        
        if ( length(findstr(tline, 'SEGMENT')) ) > 0
            SEG(s_, :) = textscan(tline, '%f %*d %*s %d %d %d %d', 1);            
            s_ = s_ + 1;
        end
    end

    fclose(fid);
    clear i b_ n_ c_ s_;
    
    % z (altitude)
    z = double(cell2mat(GPS(:, 7)));
    
    % TODO: this thing of taking 5th element is just a temporary fix
    % it should be the items which have the same time between NAV and GPS
    p = [double(cell2mat(NAV(:, 1))) ... % time
        double(cell2mat(NAV(:, 4))) ... % x ~ double(cell2mat(GPS(:, 3)))
        double(cell2mat(NAV(:, 5))) ... % y ~ double(cell2mat(GPS(:, 4)))
        z(1:5:end)];
    
    clear z;
    
    circles = [double(cell2mat(CIR(:, 1))) ... % time
               double(cell2mat(CIR(:, 2))) ... % circle center x
               double(cell2mat(CIR(:, 3))) ... % circle center y
               double(cell2mat(CIR(:, 4)))]; % circle radius
    
    segments = [double(cell2mat(SEG(:, 1))) ... % time
                double(cell2mat(SEG(:, 2))) ... % line point 1 x
                double(cell2mat(SEG(:, 3))) ... % line point 1 y
                double(cell2mat(SEG(:, 4))) ... % line point 2 x
                double(cell2mat(SEG(:, 5)))]; % line point 2 y
    
    % power
    w = double(cell2mat(BAT(:, 3)));
    w = w * (ah / 10.0);
    
    w = [double(cell2mat(BAT(:, 1)))  w];
end

