%% compute a spectrum of the pressure data you used in HW3, 4

clear all; close all;

% create empty arrays to hold time and temp data
time = [];
pressure = [];

time = [time; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-2015.nc'),'time')];
pressure = [pressure; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-2015.nc'),'pressure')];

% remove bad data using the flagged data from .nc file
pressure_flagPrimary = [];
pressure_flagPrimary = [pressure_flagPrimary; ncread(strcat('http://sccoos.org/thredds/dodsC/autoss/scripps_pier-2015.nc'),'pressure_flagPrimary')];

% looping through to remove bad data from pressure record
for i = 1:length(pressure)
    if pressure_flagPrimary(i) ~= 1
        pressure(i) = nan;
    end
end

% consider a period with equal increments

% this for the first 34 days of the 2015 record
% divide into 2 segments because need at least 14 days to resolve 12-hr tidal
% timescale

% time differences
time = double(time);
t_diff = diff(time);

% finding segment of data with even spacing
cutoff = find(t_diff(1:82236)>t_diff(1),1);
pressure_sub = pressure(1:cutoff-1); % subsampled pressure
pressure_sub = pressure_sub(1:length(pressure_sub)-3);
time_sub = time(1:cutoff-1); % subsampled times
date0=datenum(1970,1,1); % give reference date (first date)
time_sub = double(time_sub)/86400+date0; % in units of days (conversion: seconds/day)

M_pressure = 2; % number of segments
N_pressure = length(pressure_sub)/M_pressure; % datapoints per segment

% split into two segments
pressure_sub1 = reshape(pressure_sub,N_pressure,M_pressure);
pressure_sub2 = pressure_sub(N_pressure/2+1:length(pressure_sub)-N_pressure/2);
pressure_sub2 = reshape(pressure_sub2,N_pressure,M_pressure-1);
pressure_sub3 = [pressure_sub1 pressure_sub2]; % 3 segments, overlapping

pressure_sub3 = detrend(pressure_sub3).*(hann(4154)*ones(1,3));
pressure_sub3 = fft(pressure_sub3);
pressure_sub3amp = (abs(pressure_sub3(1:N_pressure/2+1,:)).^2)/N_pressure; % amplitude of first half
pressure_sub3mean = mean(pressure_sub3amp,2);
pressure_sub3mean = pressure_sub3mean'; % turn into row vector

frequency = (0:2078-1)/(2078*361)*(24*3600); 

nu = 2*3; % DOF = 2*number of segments, which seems to be best fit
err_high = nu/chi2inv(0.05/2,nu);
err_low = nu/chi2inv(1-0.05/2,nu);
ratio_chi2 = err_high/err_low;


figure
%semilogy(1:p/2-1, a_amp_mean, '-b', [p/4 p/4],[err_low err_high]*a_amp_mean(p/4), '-r')
semilogy(frequency,pressure_sub3mean, '-b', [20 20],[err_low err_high]*pressure_sub3mean(1000), '-r');
xlabel('\fontsize{14}cycles per day')
ylabel('\fontsize{14}dbar^{2}/cpd')
title('\fontsize{16}Spectrum of 2015 Scripps Pier Pressure');
legend('\fontsize{12}Pier Pressure');
