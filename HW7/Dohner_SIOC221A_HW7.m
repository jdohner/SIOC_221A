% file SIOC 221A HW 7
% 
% author Julia Dohner
%
% due date November 21, 2017
%
% I certify that this represents my own work and that I have not worked
% with classmates or other individuals to compelte this assignment. -JLD

clear all; close all;

% TODO: rename air temp and sea temp as T0 and Ta

%% inspect the data

% create empty arrays to hold time and temp data
time = [];
swTemp = [];
airTemp = [];

% time is 04/21/2015 to 06/22/2016
time = [time; ncread(strcat('OS_Stratus_2015_D_M.nc'),'TIME')];
swTemp = [swTemp; ncread(strcat('OS_Stratus_2015_D_M.nc'),'TEMP')];
airTemp = [airTemp; ncread(strcat('OS_Stratus_2015_D_M.nc'),'AIRT')];

% examining the time increments between adjacent measurements
t_diff = diff(time);
minDiff = min(t_diff); 
maxDiff = max(t_diff);
plot(t_diff);
% after plotting the values, I see that the t_diff values range between
% 0.0006944444 and 0.00069444452. The differences are sufficiently small
% to proceed with the Fourier transform.

% checking for NaN's in data:
numNaNsw = sum(isnan(swTemp(:)));
numNaNair = sum(isnan(airTemp(:)));

% plot the time series
date0=datenum(1950,1,1); % give reference date (first date)
time2 = double(time)+date0; 
figure('name','Seawater and Air Temperatures');
plot(time2, swTemp, '-b','LineWidth',1);
hold on
plot(time2, airTemp, '-r','LineWidth',1);

set(gca,'FontSize',16);
title('Stratus Seawater and Air Temperatures 04/21/2015 to 06/22/2016');
xlabel('Date');
datetick('x','mm/dd/yy')
ylabel('^{\circ}C');
legend('\fontsize{12}Sea', '\fontsize{12}Air');

%% compute the spectra for T0 and Ta

% lecture 7 is good for this

% data is taken every minute (86400 (seconds in a day)*t_diff = ~60 secs)
% split into 2 chunks to resolve seasonal cycle
M = 2; % number of segments splitting data into
N = 2*floor(length(time)/4); % length of each chunk of data (aka segment 
% length,  making sure divisible by 4)

% split into segments
swTemp2 = reshape(swTemp(1:N*2),N,M);
airTemp2 = reshape(airTemp(1:N*2),N,M);
% don't need to do overlapping segments

% apply Hanning window
swTemp2 = detrend(swTemp2).*(hann(length(swTemp2(:,1)))*ones(1,M));
airTemp2 = detrend(airTemp2).*(hann(length(airTemp2(:,1)))*ones(1,M));
T = length(airTemp2)/1440; % total time in days in each segment (1440 mins/day)
normalizationFactor = T/(N^2); 
swTemp2_amp = 2.*(abs(swTemp2(1:swTemp2/2+1,:)).^2).*normalizationFactor; % amplitude of first half

%above was done based off HW6, but I haven't gotten corrections back for
%that yet. Instead do based off Sarah Gille's HW5 solutions. 
