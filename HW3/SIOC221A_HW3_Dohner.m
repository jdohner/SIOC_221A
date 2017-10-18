% file SIOC 221A HW 3
% 
% author Julia Dohner
%
% due date October 19, 2017

clear all; close all;

numYears = 2017-2005 + 1;

%% plotting the Scripps Pier 2015 pressure record

% create empty arrays to hold time and temp data
time = [];
pressure = [];

time = [time; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-2015.nc'),'time')];
pressure = [pressure; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-2015.nc'),'pressure')];


% remove bad data using the flagged data from .nc file
pressure_flagPrimary = [];
pressure_flagPrimary = [pressure_flagPrimary; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-2015.nc'),'pressure_flagPrimary')];

%looping through to remove bad data from pressure record
for i = 1:length(pressure)
    if pressure_flagPrimary(i) ~= 1
        pressure(i) = nan;
    end
end

% examining the time increments between adjacent measurements
X = diff(time);
%figure
%plot(X);
%plot(X(1:4000))

% plot the time series
date0=datenum(1970,1,1); % give reference date (first date)
time2 = double(time)/24/3600+date0; % divide the time by 24*3600 to convert seconds into days since 1970
figure('name','Scripps_Pier_Pressure_2015');
plot(time2, pressure,'LineWidth',1);

% label the x-axis in months
set(gca,'FontSize',16);
title('Scripps Pier Pressure');
xlabel('Date');
datetick('x','mm/dd/yy')
ylabel('dbar');


% The data are for the most part uniformly spaced by 361 seconds, aside
% from later in the record where there seems to be a huge break of 17 days
% (max(X) = 1476982 seconds), and some places where the data is taken at
% closer intervals (min(X) = 220 seconds)

%% plotting just the first month of 2015

% 2592000 seconds in 30 days
% If a measurement is taken every 361 seconds, then 30 days into the record
% should be roughly the first 7180 measurements (2592000/361) in the 2015 
% series. My record is 30 days long.
time3 = time(1:7180);
pressure2 = pressure(1:7180);

% plot the time series
date0=datenum(1970,1,1); % give reference date (first date)
time4 = double(time3)/24/3600+date0; % divide the time by 24*3600 to convert seconds into days since 1970
figure('name','Scripps_Pier_Pressure_January_2015');
plot(time4, pressure2,'LineWidth',1);
% TODO: format x axis date ticks
title('Scripps Pier Pressure in Janury 2015');
xlabel('Date');
set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x','mm/dd/yy')
ylabel('dbar');

%% Least squares fit

% defining sine and cosine components of major tidal constituents

% convert period to days (to match x axis time units)
O1_sin = sin(2*pi*time4/(25.83/24)); %O1: principal lunar diurnal
O1_cos = cos(2*pi*time4/(25.83/24)); 
K1_sin = sin(2*pi*time4/(23.93/24)); %K1: luni-solar diurnal
K1_cos = cos(2*pi*time4/(23.93/24));
M2_sin = sin(2*pi*time4/(12.42/24)); %M2: principal lunar
M2_cos = cos(2*pi*time4/(12.42/24));


% TODO: do I want the time4(:) in here? I took it out of second column of A
% matrix. Got rid of the error when I removed it

A2=[ones(length(time4),1) O1_sin O1_cos K1_sin K1_cos M2_sin M2_cos];
x2=inv(A2'*A2)*A2'*pressure2;
figure('name','Pier_Pressure_Tidal_LSF');
plot(time4,A2*x2,'m','LineWidth',2)
hold on
plot(time4, pressure2,'LineWidth',1);
set(gca,'FontSize',16);
title('Pier Pressure Least Squares Fit of 3 Major Tidal Constituents');
xlabel('Time');
datetick('x','mm/dd/yy')
ylabel('TODO');


% The mean is 3.4898 (first row in x2 vector)
% 
% Total amplitude = square root of the sum of the squares of the sine and 
% cosine amplitudes)
% Units of mean and amplitude are decibars
amplitude_O1_jan = sqrt((x2(2,1))^2 + (x2(3,1))^2);
amplitude_K1_jan = sqrt((x2(4,1))^2 + (x2(5,1))^2);
amplitude_M2_jan = sqrt((x2(6,1))^2 + (x2(7,1))^2);

%% Stationarity of the tide

% repeating the least squares fit for 30 days roughly near August 2015 

% Starting 7/12 of the way through the time record (82237 measurements).
% If a measurement is taken every 361 seconds, then 30 days into the record
% should be roughly the next 7180 measurements (2592000/361). My record is 
% 30 days long.
kAugustStart = floor((7/12)*82237); 
kAugustEnd = kAugustStart + 7180; 
timeAugust = time(kAugustStart:kAugustEnd, 1);
timeAugust = double(timeAugust)/24/3600+date0; % in units of days
pressureAugust = pressure(kAugustStart:kAugustEnd, 1);

% defining sine and cosine components of major tidal constituents

% convert period to days (to match x axis time units)
O1_sin_aug = sin(2*pi*timeAugust/(25.83/24)); %O1: principal lunar diurnal
O1_cos_aug = cos(2*pi*timeAugust/(25.83/24)); 
K1_sin_aug = sin(2*pi*timeAugust/(23.93/24)); %K1: luni-solar diurnal
K1_cos_aug = cos(2*pi*timeAugust/(23.93/24));
M2_sin_aug = sin(2*pi*timeAugust/(12.42/24)); %M2: principal lunar
M2_cos_aug = cos(2*pi*timeAugust/(12.42/24));


A3=[ones(length(timeAugust),1) O1_sin_aug O1_cos_aug K1_sin_aug K1_cos_aug M2_sin_aug M2_cos_aug];
x3=inv(A3'*A3)*A3'*pressureAugust;
figure('name','Pier_Pressure_Tidal_LSF_Summer');
plot(timeAugust,A3*x3,'m','LineWidth',2)
hold on
plot(timeAugust, pressureAugust,'LineWidth',1);
set(gca,'FontSize',16);
title('Pier Pressure Least Squares Fit of 3 Major Tidal Constituents - Summer');
xlabel('Time');
datetick('x','mm/dd/yy')
ylabel('TODO');

% The mean is 3.4898 (first row in x2 vector)
% 
% Total amplitude = square root of the sum of the squares of the sine and 
% cosine amplitudes)
% Units of mean and amplitude are decibars
amplitude_O1_aug = sqrt((x3(2,1))^2 + (x3(3,1))^2);
amplitude_K1_aug = sqrt((x3(4,1))^2 + (x3(5,1))^2);
amplitude_M2_aug = sqrt((x3(6,1))^2 + (x3(7,1))^2);

% The tidal amplitudes for the K1 tide decreases by 0.07. The amplitude of 
% M2 decreases by 0.01, and the amplitude of O1 does not change. The tidal 
% component that changes the most is the K1 tide. Because the earth rotates 
% on a tilted axis relative to our axis of rotation around the sun, the
% distance of the pier from the sun changes on a seasonal scale, changing
% the influence of solar forcing on local tides. The K1 tide is the only of
% the three tides that contains a component reflecting the influence of the
% sun on tides, thus it makes sense that K1 shows the biggest change in
% amplitude between January and August.

%% Chi squared and the misfit. (Good name for a short story)

% y is vector containing pressure data (pressureAugust)
% A is the matrix containing ones, sines and cosines (A3)
% x is matrix containing mean and amplitudes (x3)
% TODO: Question: is sigma the uncertainty for each point? Or sigma for the
% whole data?
sigma = std(pressureAugust);
%sum = symsum((((pressureAugust - A3*x3).^2)/sigma.^2), 1, length(pressureAugust));

chiSquared = 0; % initialized chiSquared variable before loop
for k = 1:length(pressureAugust)
    chiSquared = chiSquared + (pressureAugust(k) - (A3*x3)









