
% load data from pprz log
[p w segments circles] = loadpprz('18_12_12__09_53_01_SD.data', 2.5);

% filter data about the cruise (~490.5->1226)
p(p(:, 1) < 680, :) = [];
p(p(:, 1) > 1000, :) = [];

w(w(:, 1) < 680, :) = [];
w(w(:, 1) > 1000, :) = [];

segments(segments(:, 1) < 680, :) = [];
segments(segments(:, 1) > 1000, :) = [];

circles(circles(:, 1) < 680, :) = [];
circles(circles(:, 1) > 1000, :) = [];

% normalizing and merging
% units are now: x, y (utm normalized in cm), altitude (cm)
% so it's expressed in cm, meaning this is the way to get in in meters
% data from NAV, so no more normalization needed... Except for height
%p(:, 2) = p(:, 2) - min(p(:, 2));
%p(:, 3) = p(:, 3) - min(p(:, 3));
%p(:, 2) = p(:, 2) / 100;
%p(:, 3) = p(:, 3) / 100;
p(:, 4) = p(:, 4) / 100;

%% animation
delay = p(2, 1) - p(1, 1);
    
figure(1);
 
for i=1:size(p(:, 1))
    
    subplot(2, 2, [1 2]);
    plot3(p(1:i, 2), p(1:i, 3), p(1:i, 4))
    %title(datestr(p(i, 1) / (24 * 60 * 60), 'HH:MM:SS'));
    
    xlim([min(p(:, 2)) max(p(:, 2))]);
    ylim([min(p(:, 3)) max(p(:, 3))]);
    zlim([0 max(p(:, 4))]);
    
    subplot(2, 2, 3);
    plot(p(1:i, 2), p(1:i, 3))
    %title(p(i, 1));
    
    xlim([min(p(:, 2)) max(p(:, 2))]);
    ylim([min(p(:, 3)) max(p(:, 3))]);
    
    if any(fix(circles(:, 1)) == fix(p(i, 1)))
        subplot(2, 2, 3);
        syms x y;
        circle = circles(fix(circles(:, 1)) == fix(p(i, 1)), 2:end);
        
        f = (x - circle(1))^2 + (y - circle(2))^2 - circle(3)^2;
        
        subplot(2, 2, 4);
    
        fimplicit(f);
        xlim([min(p(:, 2)) max(p(:, 2))]);
        ylim([min(p(:, 3)) max(p(:, 3))]);
    end
    
    if any(fix(segments(:, 1)) == fix(p(i, 1)))
        subplot(2, 2, 3);
        syms x y;
        segment = segments(fix(segments(:, 1)) == fix(p(i, 1)), 2:end);
        
        if segment(1) == segment(3)
            f = x - segment(1);
        else
            f = (y - segment(2)) / (x - segment(1)) - (segment(4) - segment(2)) / (segment(3) - segment(1));
        end
        
        subplot(2, 2, 4);
    
        fimplicit(f, [min(p(:, 2)) max(p(:, 2)) min(p(:, 3)) max(p(:, 3))]);
    end
    
    pause(delay)
    
end
 
 