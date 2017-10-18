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
figure
plot(X(1:4000))

% plot the time series
date0=datenum(1970,1,1); % give reference date (first date)
time2 = double(time/24/3600+date0); % divide the time by 24*3600 to convert seconds into days since 1970
figure('name','Scripps_Pier_Pressure_2015');
plot(time2, pressure,'LineWidth',1);

% label the x-axis in months
datetick('x','yyyy');
set(gca,'FontSize',16);
title('Scripps Pier Pressure');
xlabel('Date');
datetick('x','mm/dd/yy','keeplimits')
ylabel('dbar');


% The data are for the most part uniformly spaced by 361 seconds, aside
% from later in the record where there seems to be a huge break of 17 days
% (max(X) = 1476982 seconds), and some places where the data is taken at
% closer intervals (min(X) = 220 seconds)

%% plotting just the first month of 2015

% 2592000 seconds in 30 days
firstTime = min(time);
% if measurement is taken every 361 seconds, then 30 days into the record
% should be roughly the first 7180 measurements (2592000/361) in the 2015 
% series 
time3 = time(1:7180);
pressure2 = pressure(1:7180);

% plot the time series
date0=datenum(1970,1,1); % give reference date (first date)
time4 = double(time3/24/3600+date0); % divide the time by 24*3600 to convert seconds into days since 1970
figure('name','Scripps_Pier_Pressure_January_2015');
plot(time4, pressure2,'.');

% label the x-axis in months
datetick('x','yyyy');
set(gca,'FontSize',16);
title('Scripps Pier Pressure in Janury 2015');
xlabel('Date');
datetick('x','mm/dd/yy')
ylabel('dbar');




