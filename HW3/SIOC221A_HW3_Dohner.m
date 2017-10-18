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
firstTime = min(time);
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
% label the x-axis in months
%datetick('x','yyyy');
%set(gca,'FontSize',16);
%xlinspace = linspace(min(time4), max(time4), 5);
title('Scripps Pier Pressure in Janury 2015');
xlabel('Date');
set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x','mm/dd/yy')
ylabel('dbar');

%% Least squares fit
%TODO: are we only supposed to fit our one month of data for this? or
%entire year?

% % assuming time matches the SST standard with days counting from
% % January 1, 1800
% % and data is the monthly sea surface temperature
% plot(time4,pressure2,'LineWidth',2); hold on
% A=[ones(length(time),1) time(:)];
% x=inv(A'*A)*A'*data;
% plot(time+datenum(1800,1,1),A*x,'r','LineWidth',2)

%time5 = double(time3)/24/3600+date0;

% defining sine and cosine components of major tidal constituents

% need the units of the period to be the same as the x axis time units
O1_sin = sin(2*pi*time4/25.83); %O1: principal lunar diurnal
O1_cos = cos(2*pi*time4/25.83); 
K1_sin = sin(2*pi*time4/23.93); %K1: luni-solar diurnal
K1_cos = cos(2*pi*time4/23.93);
M2_sin = sin(2*pi*time4/12.42); %M2: principal lunar
M2_cos = cos(2*pi*time4/12.42);

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

% LSF should basically match the data
% whats the amplitude of the tide components? can calculate from sin & cos
% mean ~3.5
% tidal amplitudes will be 1-2

%TODO: what is the mean and what are the total amplitudes of the three
%tidal constituents? (Total amplitude should be determined from the square
%root of the sum of the squares of the sine and cosine amplitudes)

