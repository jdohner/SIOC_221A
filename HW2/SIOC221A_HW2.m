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

% retrieve flagged data
for i = 1:numYears
    n = (i-1) + 2005; % get value for each year starting with first year
    salinity_flagPrimary = [salinity_flagPrimary; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'salinity_flagPrimary')];

end


%looping through to remove bad data from salinity record
for i = 1:length(salinityExtend)
    if salinity_flagPrimary(i) ~= 1
        salinityExtend(i) = nan;
    end
end


% plot the time series
date0=datenum(1970,1,1); % give reference date (first date)
time = double(timeExtend/24/3600+date0);
figure('name','Scripps_Pier_Salinity_2005-2017_Automated');
plot(time, salinityExtend,'LineWidth',1);

% label the x-axis in months
datetick('x','yyyy');
set(gca,'FontSize',16);
title('Scripps Pier Salinity - Automated Record');
xlabel('Date');
datetick('x','mm/dd/yy','keeplimits')
ylabel('Practical Salinity Units (1e-3)');

% mean salinity
meanSalinityAuto = nanmean(salinityExtend);
stdSalinityAuto = nanstd(salinityExtend);

%% plotting the manual salinity record

filein = 'SIO_SALT_1916-201410.txt';
fileID = fopen(filein);
%read in header info
headerInfo = textscan(fileID, '%s', 27, 'delimiter', '\n');
dataInfo = textscan(fileID, '%s', 9, 'delimiter', '\t');
% read in the data
C = textscan(fileID, '%f %f %f %f %f %f %f %f %f', 'delimiter', '\t \t \t \t \t \t \t \t \t');

%extract the data from the cell matrix
yearData = C{1};
monthData = C{2};
dayData = C{3};
timeData = C{4}; % seems like these are all NaN's
timeFlagData = C{5};
salinityData = C{6};
salinityFlagData = C{7};
bottleSalData = C{8};
bottleSalFlagData = C{9};


%looping through to remove bad data from salinity record
for i = 1:length(yearData)
% ignoring flags for time
%     if timeFlagData(i) ~= 0
%         salinityData(i) = nan;
    if salinityFlagData(i) ~= 0
            salinityData(i) = nan;
    elseif bottleSalFlagData(i) ~= 0
        bottleSalData(i) = nan;
    end
end

% convert all nan times to noon
for i = 1:length(timeData)
    if isnan(timeData(i)) == 1
        timeData(i) = 1200;
    end
end

% convert time data to hours and minutes
time_hour = floor(timeData/100);
time_minute = timeData-time_hour*100;
time_second = zeros(length(yearData),1);        
        
%turn the year, mo, day, time into a MATLAB date
dateData = datenum(yearData, monthData, dayData, time_hour, time_minute, time_second);

%plot time series
figure('name','Scripps_Pier_Salinity_1916-2004_Manual');
plot(dateData,salinityData,'-')
set(gca,'FontSize',16);
t1 = datenum('22-august-1916');
t2 = datenum('31-october-2014');
xlim([t1 t2]);
%label the plot
datetick('x','mm/dd/yy','keeplimits')
xlabel('Date')
title('Scripps Pier Salinity - Manual Record')
ylabel('Practical Salinity Units s(1e-3)');

% mean salinity
meanSalinityManual = nanmean(salinityData);
stdSalinityManual = nanstd(salinityData);

%% subsampling to compare manual and auto

% only compare data from manual and auto that were taken at the same time
% comparing between Jan 1 2005 and Jan 1 2015


% subAuto = nan[];
% % loop from Jan 1 2005 to Jan 1 2015
% startDate = datenum('2005');
% endDate = datenum('2015');
% kStartManual = find(yearData == startDate);
% for i = startDate:endDate

testDateVector = datenum(yearData, monthData, dayData);
startDate = ('1-january-2005');
endDate = ('1-january-2015');
kStartManual = find(testDateVector == startDate);

%% theoretical PDFs

x = meanSalinityAuto-(4*stdSalinityAuto):0.01:meanSalinityAuto+(4*stdSalinityAuto);

% automated data pdf for a variety of distributions (given the mean & std)
%pd1 = makedist('Normal', meanSalinityAuto, stdSalinityAuto);

%y = pdf(pd1,x);

%y = normpdf(x,meanSalinityAuto,stdSalinityAuto);


gaussianY = pdf('Normal', x, meanSalinityAuto, stdSalinityAuto);

%std = square root of variance
% variance = std^2
uniformY = pdf('Uniform', x, meanSalinityAuto, stdSalinityAuto);
bimodalY = pdf('Binomial', x, meanSalinityAuto, stdSalinityAuto);

figure
subplot(3,1,1);
plot(x,gaussianY);
subplot(3,1,2);
plot(x,uniformY);
subplot(3,1,3);
plot(x,bimodalY);



