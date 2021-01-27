
% simulation from the equations from Hector

delta = 0.01;
time = delta;

m = 1; % mass of the aircraft [kg]
W = m*9.8; % weight force [N]

cl = 9.8/15^2; % coefficient to be determined experimentally
cth = 15/50; % same

kp = 1; % positive gain constant to be also determined experimentally

th = 50; % nominal vlalue of the throttle

hd = 25; % desired altitude [m]

w = [1 0]; % wind vector (x, y axis velocity

vv = [15]; % vertical velocity vector (initial 15 [m/s])
h = [25]; % altitude vector (inital 25 [m])

p = [0; 0]; % position in space (initial 0,0 [m])
theta = [0.1]; % initial angle

store = []; % just storing all the records for debugging

while true
    
    dth = kp*(hd-h(end)); % dth is the change of throttle we use to adjust 
                          % the altitude
    
    vs = cth*(th+dth)-w*[cos(theta(end)); sin(theta(end))]; % airspeed of  
                                                            % the aircraft 
                                                            % [m/s]
    % Hector original formula was wbx; I changes it in -w*[cos(theta);
    % sin(theta)] which substract the contribution of the wind
    
    L = cl*vs^2;
    
    av = (L-W)/m; % vertical accelaration [m/s^2]
    
    % verical dynamics
    h = [h; h(end)+vv(end)*delta]; % integration with just Euler
    vv = [vv; vs+av*delta];
    
    % you get the pdangle here in the old simulation instead of theta...
    u = theta(end);
    
    % horizontal kinematics
    p = [p p(:,end)+(cth*(th+dth)*...
           [cos(theta(end)); sin(theta(end))]+w)*delta];
    theta = [theta; theta(end)+u*delta];
    
    
    store = [store; dth vs L av u];
    
    time = time + delta;
    
    if time >= 60 % 1 minute? It's over...
        break;
    end
    
end

