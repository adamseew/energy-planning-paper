

% load data from simulation 3 D (raw5)

%load('../../data/simulation3/raw5/updated/simulation3D.mat');

% load data from simulation 

load('../../data/simulation3/raw2/updated/simulation3A.mat');

qq0 = q(:, 200/delta);
yy = [];

ddelta = delta/5; % refinining the delta little more

AAd = A*ddelta+eye(2*r+1);

% get the estimated coefficients with KF at time 200 s on
% that's the horizon of the plot in the results section in the papaer

yyy = [];
yyr = [];

jj = 1;

for j=20000:period/delta:length(q)
    
    qq0 = q(:, uint64(j));
    yy = [];
    
    for ii=0:ddelta:200 % evolution up to the next 200 s

        qq1 = AAd*qq0;
        yy = [yy; C*qq1];
        qq0 = qq1;
    
    end
    
    yyy = [yyy yy];
    yyr = [yyr imresize(yy, [400 1])];
    
    figure(1)
	plot(linspace(200, 400, length(yy)), yy)
   
    legendInfo{jj} = [num2str(uint64(j))]; 
    
    hold on;

    figure(2)
	plot(linspace(200, 400, length(yyr(:,jj))), yyr(:,jj))
    
    legendInfo{jj} = [num2str(uint64(j))]; 
    jj = jj + 1;
    
    hold on;
    
end

legend(legendInfo);

csvwrite('evolution_simulation3A.csv', [linspace(200, 400, 400)' yyr]);

