
% matlab script to generate the plots in the illustrative abstract

% read original data
x = csvread('../../data/simulation1/pprz_throttle.csv');
x = x(:,3); % the energy is in the 3rd column
x2 = resample(x, 100, 640); % shrink data
x2(1) = x(1); % make sure first value is the same in both

% the plot with the energy it time
t = linspace(0, 320, 100)';
figure(1);
plot(t, x2)
xlabel('Time (sec)');
ylabel('Power (W)');

% spectral analysis
period = 16; % this value is observed from data
             % I looked where are the peaks
xline(50);
xline(100);
xline(150);
xline(200);
xline(250);
xline(300); % just some lines at the hypothetical period 

fs = 1/period; % frequency
y = fft(x - mean(x)); % discrete fourier transform of the energy signal

n = length(x); % number of samples
f = (0:n-1)*(fs/n); % frequency range
power = abs(y).^2/n; % power of the DFT

figure(2);
plot(f,power)
xlabel('Frequency')
ylabel('Power')

% power spectrum centered at 0 frequency (more insightful)
y0 = fftshift(y);         % shift y values
power0 = abs(y0).^2/n;    % 0-centered power
n = 100; % taking just 10 frequencies
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range

figure(3);
plot(f0,power0((end-n)/2+1:(end+n)/2))
xlabel('Frequency')
ylabel('Power')

csvwrite('periodic_energy.csv', [t x2]);
csvwrite('spectrum.csv', [f0' power0((end-n)/2+1:(end+n)/2)]);

