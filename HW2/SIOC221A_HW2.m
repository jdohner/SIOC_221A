% file SIOC 221A HW 2
% 
% author Julia Dohner
%
% due date October 12, 2017
%
% file Plotting salinity record from Shores Station Program at Scripps Pier

clear all; close all;

numYears = 2017-2005 + 1;

%% plotting the automated salinity record

% create empty arrays to hold time and temp data
timeExtend = [];
salinityExtend = [];

% loop through each year starting with 2005 to retrieve data
for i = 1:numYears
    n = (i-1) + 2005; % get value for each year starting with first year
    timeExtend = [timeExtend; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'time')];
    salinityExtend = [salinityExtend; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'salinity')];

end

% remove bad data using the flagged data from .nc file
% create empty arrays to hold time and temp data
salinity_flagPrimary = [];
salinity_flagSecondary = [];

for i = 1:numYears
    n = (i-1) + 2005; % get value for each year starting with first year
    salinity_flagPrimary = [salinity_flagPrimary; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'salinity_flagPrimary')];
    salintiy_flagSecondary = [salinity_flagSecondary; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'salinity_flagSecondary')];

end

% %TODO: find indices of bad data points

%looping through to remove bad data from salinity recor
for i = 1:length(salinityExtend)
    if salinity_flagPrimary(i) ~= 1
        salinityExtend(i) = nan;
%     TODO: not yet sure what to do with salinity_flagSecondary
%     else if salintiy_flagSecondary(i) ~= 1
%             salinityExtend(i) = nan;
%     end
    end
end


% plot the time series
date0=datenum(1970,1,1); % give reference date (first date)
time = double(timeExtend/24/3600+date0);
figure('name','Scripps_Pier_Salinity_2005-2017_Automated');
plot(time, salinityExtend,'LineWidth',1);

% % plot the flagged values
% date0=datenum(1970,1,1); % give reference date (first date)
% time = double(timeExtend/24/3600+date0);
% figure('name','Scripps_Pier_Salinity_2005-2017_Automated');
% plot(time, salinity_flagPrimary,'LineWidth',1);

% label the x-axis in months
datetick('x','yyyy');
set(gca,'FontSize',16);
title('Scripps Pier Salinity - Automated Record');
xlabel('year','FontSize',16);
ylabel('Salinity (1e-3)', 'FontSize',16);

%% plotting the manual salinity record
