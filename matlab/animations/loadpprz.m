
function [p, w, segments, circles] = loadpprz(fn, ah)

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
            GPS(i, :) = textscan(tline, '%f %*d %*s %d %d %d %d %d %d %d %d %d %d', 1);
            i = i + 1;
        end
        
        if  (length(findstr(tline, 'NAVIGATION_REF'))) > 0
            continue;
        end
        
        if ( length(findstr(tline, 'NAVIGATION')) ) > 0
            NAV(n_, :) = textscan(tline, '%f %*d %*s %d %d %d %d %d %d %d %d', 1);
            n_ = n_ + 1;
        end
   
        if ( length(findstr(tline, 'BAT')) ) > 0
            BAT(b_, :) = textscan(tline, '%f %*d %*s %d %d %d %d %d %d %d %d', 1);            
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
    clear i _b;
    
    % x
    % x = double(cell2mat(GPS(:, 3)));
    % or
    x = double(cell2mat(NAV(:, 4)));
    % y
    % y = double(cell2mat(GPS(:, 4)));
    % or
    y = double(cell2mat(NAV(:, 5)));
    % z (altitude)
    z = double(cell2mat(GPS(:, 7)));
    
    % TODO this thing of taking 5th element is just a temporary fix
    % it should be the items which have the same time between NAV and GPS
    p = [(double(cell2mat(NAV(:, 1))) - double(cell2mat(NAV(1, 1)))) x y z(1:5:end)];
    
    % circle center E
    crc_e = double(cell2mat(CIR(:, 2)));
    % circle center N
    crc_n = double(cell2mat(CIR(:, 3)));
    % circle radius
    crc_rad = double(cell2mat(CIR(:, 4)));
    
    circles = [(double(cell2mat(CIR(:, 1))) - double(cell2mat(CIR(1, 1)))) crc_e crc_n crc_rad];
        
    % segment E
    sgm_e1 = double(cell2mat(SEG(:, 2)));
    % segment N (first point)
    sgm_n1 = double(cell2mat(SEG(:, 3)));
    % segment E
    sgm_e2 = double(cell2mat(SEG(:, 4)));
    % segment N (second point)
    sgm_n2 = double(cell2mat(SEG(:, 5)));
    
    segments = [(double(cell2mat(SEG(:, 1))) - double(cell2mat(SEG(1, 1)))) sgm_e1 sgm_n1 sgm_e2 sgm_n2];
    
    % power [w]
    w = double(cell2mat(BAT(:, 3)));
    w = w * (ah / 10.0);
    
    w = [(double(cell2mat(BAT(:, 1))) - double(cell2mat(BAT(1, 1)))) w];
end

