% file SIOC 221A HW 1
% 
% author Julia Dohner, borrowed heavily from Sarah Gille
%
% date October 2nd, 2017

clear all; close all; clc

% access sccoos data via THREDDS
ncid=netcdf.open('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-2017.nc','NOWRITE');

% determine the variable number that you want to read
ntime = netcdf.inqVarID(ncid,'time')
ntemp = netcdf.inqVarID(ncid,'temperature')
% these will return 0 and 1, telling you that time is variable number 0 and temperature is variable number 1

% read data
pierTime = netcdf.getVar(ncid,ntime);
pierSST = netcdf.getVar(ncid,ntemp);

% note that time is measured in seconds since January 1, 1970
% so define a reference date
date0=datenum(1970,1,1);

%% plot the time series

% use double to force the time to be a real number
% divide the time by 24*3600 to convert seconds into days since 1970
time = double(pierTime/24/3600+date0);
figure('name','Scripps_Pier_SST_2017');
plot(time,pierSST,'LineWidth',1)

% label the x-axis in months
datetick('x','mmm');
set(gca,'FontSize',16);
title('Scripps Pier SST in 2017');
xlabel('months (of 2017)','FontSize',16);
ylabel('temperature (?oC)', 'FontSize',16);

% From the plot there appears to be a steady upward trend of SST , notably
% from March until the present, with a large peak in values at the end of
% July. There also seem to be some temperature swings in the months prior.
% The temperature appears to be steadiest between January and March.

%% compute the mean and standard deviation

meanSST = mean(pierSST);
stdSST = std(pierSST);

% these statistics tell us what value we should expect if we took the SST 
% many times at the pier (mean) and whether a deviation from this value is
% typical (indicated by whether a measurement falls within one (or however
% many) standard deviations from the mean). They can tell us what the
% average SST was in 2017, which may be useful when compared to the average
% SST from other years at the pier. 

%% compute an empirical probability density function for SST

figure('name','PDF_Scripps_Pier_SST_2017');
histogram(pierSST,'Normalization','pdf');

% label the x-axis in temperature
set(gca,'FontSize',16);
title('Probability Density Function of SST');
xlabel('temperature (oC)','FontSize',16);
ylabel('probability', 'FontSize',16);

% The plot does not look like any of the distributions we discussed in
% class.

%% extending the record
% extending record to include SST data from 2005 to 2017



% access sccoos data via THREDDS

numYears = 2017-2005 + 1;
ncidArray = zeros(1,numYears); % create empty array to hold ncids for each file

for i = 2005:2017
    % get file names for each of timeseries to be called from sccoos
    fileName = strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-', num2str(i), '.nc');
    for j=1:numYears
        ncidArray(j) = netcdf.open(fileName,'NOWRITE');
    end
    
end

% determine the variable number that you want to read
ntimeExtend = netcdf.inqVarID(ncidArray(1),'time')
ntempExtend = netcdf.inqVarID(ncidArray(1),'temperature')
% these will return 0 and 1, telling you that time is variable number 0 and temperature is variable number 1

% empty arrays to hold time and temp data for each year between 2005-2017
pierTimeExtend = [];
pierTempExtend = [];

for i = 1:numYears;
    % get values for time and temp for element i in ncidArray
    newTime = netcdf.getVar(ncidArray(i),ntimeExtend);
    newTemp = netcdf.getVar(ncidArray(i),ntempExtend);
    % append values to time and temp arrays
    pierTimeExtend = [pierTimeExtend; newTime];
    pierTempExtend = [pierTempExtend; newTemp];
end 

% plot the time series by looping through pierTimeExtend and pierTempExtend
figure('name','Scripps_Pier_SST_2005-2017');
date0=datenum(1970,1,1); % give reference date (first date)
time = double(pierTimeExtend/24/3600+date0);
plot(time, pierTempExtend,'LineWidth',1);


% label the x-axis in months
datetick('x','mmm');
set(gca,'FontSize',16);
title('Scripps Pier SST');
xlabel('months','FontSize',16);
ylabel('temperature (oC)', 'FontSize',16);


