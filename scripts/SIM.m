%SIM
% Sim merges all the scripts and builds a complete simulation of the
% position and the energy

% data about the trajectory

answer = inputdlg(...
    {'vehicle speed [m/s]:', 'wind speed [m/s]:', ...
     'vehicle direction [angles]:', 'wind direction [angles]:', ...
     'start coordinate x [m]:', 'start coordinate y [m]:', ...
     'max power [W]:', 'min power [W]:', 'triggering point radius [m]' ...
    }, ...
    'Trajectory initialization', [1 40], ...
    {'20', '5', '270', '90', '-100', '220', '60', '30', '10'}); % asking for initial data

if isempty(answer)
    strp = [20; 5; 270; 90; -100; 220; 60; 30; 10]; % default initial data
else
    strp = str2double(answer);
end

clear answer;

% data about the energy model

answer = inputdlg(...
    {'order r:', 'characteristic time xi:', ...
     'epsilon:', 'horizon N:'
    }, ...
    'Model and algorithm initialization', [1 40], ...
    {'3', '10', '1', '10'}); % asking for initial data

if isempty(answer)
    strp2 = [3; 10; 1; 10]; % default initial data
else
    strp2 = str2double(answer);
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
si = []; % the QoS set per each i
% if there are multiple instances of si (multiple task), this operation has
% to be done again for the others
i = -1;
while true
    thisline = fgetl(fid);
    if ~ischar(thisline) 
        break; 
    elseif isempty(thisline)
        continue;
    end  %end of file
    
    if contains(thisline, "[")
        thisline(1:1) = [];
        thisline(end:end) = [];
        thisline = split(thisline, ',');
        si = [si; i, i, str2double(thisline(1)), str2double(thisline(2))];
    elseif contains(thisline, ",")
        thisline = split(thisline, ",");
        trigger = [trigger; str2double(thisline(1)) str2double(thisline(2))];
    else
        varphi = [varphi; convertCharsToStrings(thisline)];
        i = i + 1;
    end
end

%TODO: the former mechanical components of the control
ci = zeros(size(i, 2), 1);
  
fclose(fid);

% The value of the model from the modeling tool (powprof); former
% computational energy.
disp('[   $ ] Input computational energy data (from powprof) [csv]');
[bn, folder] = uigetfile('.csv');
if bn == 0
    file2 = csvread('../data/simulation3/computational_energy_2.csv');
else
    file2 = csvread(fullfile([folder, bn]));
end

clear bn folder;

addpath(genpath('energy'));

% The former computational model constraints set (per each time step all
% possibilities)
si2 = build_ck(si, file2, 0:1:i);

figure;

[k pos pow est ua] = simulate(varphi, trigger, strp, strp2, si2, ci);

