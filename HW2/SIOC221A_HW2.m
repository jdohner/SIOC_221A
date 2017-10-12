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
timeAuto = [];
salinityAuto = [];

% loop through each year starting with 2005 to retrieve data
for i = 1:numYears
    n = (i-1) + 2005; % get value for each year starting with first year
    timeAuto = [timeAuto; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'time')];
    salinityAuto = [salinityAuto; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'salinity')];

end

% remove bad data using the flagged data from .nc file
% create empty arrays to hold time and temp data
salinity_flagPrimary_Auto = [];

% retrieve flagged data
for i = 1:numYears
    n = (i-1) + 2005; % get value for each year starting with first year
    salinity_flagPrimary_Auto = [salinity_flagPrimary_Auto; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(n), '.nc'),'salinity_flagPrimary')];

end


%looping through to remove bad data from salinity record
for i = 1:length(salinityAuto)
    if salinity_flagPrimary_Auto(i) ~= 1
        salinityAuto(i) = nan;
    end
end


% plot the time series
date0=datenum(1970,1,1); % give reference date (first date)
time = double(timeAuto/24/3600+date0);
figure('name','Scripps_Pier_Salinity_2005-2017_Automated');
plot(time, salinityAuto,'LineWidth',1);

% label the x-axis in months
datetick('x','yyyy');
set(gca,'FontSize',16);
title('Scripps Pier Salinity - Automated Record');
xlabel('Date');
datetick('x','mm/dd/yy','keeplimits')
ylabel('Practical Salinity Units (1e-3)');

% mean salinity
meanSalinityAuto = nanmean(salinityAuto);
stdSalinityAuto = nanstd(salinityAuto);

%% plotting the manual salinity record

filein = 'SIO_SALT_1916-201410.txt';
fileID = fopen(filein);
%read in header info
headerInfo = textscan(fileID, '%s', 27, 'delimiter', '\n');
dataInfo = textscan(fileID, '%s', 9, 'delimiter', '\t');
% read in the data
C = textscan(fileID, '%f %f %f %f %f %f %f %f %f', 'delimiter', '\t \t \t \t \t \t \t \t \t');

%extract the data from the cell matrix
yearDataManual = C{1};
monthDataManual = C{2};
dayDataManual = C{3};
timeDataManual = C{4}; 
timeFlagDataManual = C{5};
salinityDataManual = C{6};
salinityFlagDataManual = C{7};
bottleSalDataManual = C{8};
bottleSalFlagDataManual = C{9};


%looping through to remove bad data from salinity record
for i = 1:length(yearDataManual)
% ignoring flags for time
%     if timeFlagData(i) ~= 0
%         salinityData(i) = nan;
    if salinityFlagDataManual(i) ~= 0
            salinityDataManual(i) = nan;
    elseif bottleSalFlagDataManual(i) ~= 0
        bottleSalDataManual(i) = nan;
    end
end

% convert all nan times to noon
for i = 1:length(timeDataManual)
    if isnan(timeDataManual(i)) == 1
        timeDataManual(i) = 1200;
    end
end

% convert time data to hours and minutes
time_hour = floor(timeDataManual/100);
time_minute = timeDataManual-time_hour*100;
time_second = zeros(length(yearDataManual),1);        
        
%turn the year, mo, day, time into a MATLAB date
dateData = datenum(yearDataManual, monthDataManual, dayDataManual, time_hour, time_minute, time_second);


%plot time series
figure('name','Scripps_Pier_Salinity_1916-2004_Manual');
plot(dateData,salinityDataManual,'-')
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
meanSalinityManual = nanmean(salinityDataManual);
stdSalinityManual = nanstd(salinityDataManual);

%% subsampling to compare manual and auto

% only compare data from manual and auto that were taken at the same time
% comparing between Jan 1 2005 and Jan 1 2015


% subAuto = nan[];
% % loop from Jan 1 2005 to Jan 1 2015
% startDate = datenum('2005');
% endDate = datenum('2015');
% kStartManual = find(yearData == startDate);
% for i = startDate:endDate

testDateVector = datenum(yearDataManual, monthDataManual, dayDataManual);
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

% solved for upper and lower bounds in wolfram alpha using:
% 1/12*(upper-lower)^2 = std^2
% 0.5*(upper+lower) = mean
pdUniform = makedist('Uniform','lower', 32.6138, 'upper', 34.058);
uniformY = pdf(pdUniform,x);

% solved for N and P in wolfram alpha using:
% mean = N*P
% variance = N*P(1-P)
N = 34;
P = 0.994781;
pdBimodal = makedist('Binomial',N,P);
bimodalY = pdf(pdBimodal,x);
%bimodalY = pdf('Binomial', x, meanSalinityAuto, stdSalinityAuto);

figure
subplot(3,1,1);
plot(x,gaussianY);
subplot(3,1,2);
plot(x,uniformY);
subplot(3,1,3);
plot(x,bimodalY);

%% empirical probability density functions

figure('name','PDF_Scripps_Pier_Salinity');
subplot(2,1,1)
histogram(salinityAuto,'Normalization','pdf');

set(gca,'FontSize',16);
title('Probability Density Function of SST for Automated Data');
xlabel('Salinity in PSU','FontSize',16);
ylabel('probability', 'FontSize',16);

subplot(2,1,2)
histogram(salinityDataManual,'Normalization','pdf');

minManual = nanmin(salinityDataManual);
maxManual = nanmax(salinityDataManual); %why is this 720??

set(gca,'FontSize',16);
title('Probability Density Function of SST for Manual Data');
xlim([29.64, 35]) % MATLAB wouldn't take variables here but they're the min and max of manual salinity values
xlabel('Salinity in PSU','FontSize',16);
ylabel('probability', 'FontSize',16);

