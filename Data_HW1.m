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

