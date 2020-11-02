%SIM
% Sim merges all the scripts and builds a complete simulation of the
% position and the energy

answer = inputdlg(...
    {'vehicle speed [m/s]:', 'wind speed [m/s]:', ...
     'vehicle direction [angles]:', 'wind direction [angles]:', ...
     'start coordinate x [m]:', 'start coordinate y [m]:', ...
     'max power [W]:', 'min power [W]:' ...
    }, ...
    'Initialization', [1 40], ...
    {'15', '5', '0', '0', '-100', '240', '40' '30'}); % asking for initial data

if isempty(answer)
    strp = [15; 5; 0; 0; -100; 240; 40; 30]; % default initial data
else
    strp = str2double(answer);
end

clear answer;

% retriving paparazzi log
disp('[   $ ] Input mission specification [mis]');
[bn, folder] = uigetfile('.mis');
if bn == 0
    fid = fopen('../data/simulation3/agriculture.mis', 'rt');
else
    fid = fopen(fullfile([folder, bn]), 'rt');
end

    varphi = []; % TEE per each i
trigger = []; % triggering points from i to i+1
while true
    thisline = fgetl(fid);
    if ~ischar(thisline) 
        break; 
    elseif isempty(thisline)
        continue;
    end  %end of file
    
    if contains(thisline, ",")
        thisline = split(thisline, ",");
        trigger = [trigger; str2double(thisline(1)) str2double(thisline(2))];
    else
        varphi = [varphi; convertCharsToStrings(thisline)];
    end
    
end
  
fclose(fid);

simulate(varphi, trigger, strp(1),  strp(2), strp(3), strp(4), strp(5), ...
    strp(6), strp(7), strp(8));

