% file SIOC 221A HW 4
% 
% author Julia Dohner
%
% due date October 26, 2017

clear all; close all;


%% gathering 2015 Pier pressure record

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


%% Pier pressure record for a summer month of 2015

% 30 days roughly near August 2015 

% Starting 7/12 of the way through the time record (82237 measurements).
% If a measurement is taken every 361 seconds, then 30 days into the record
% should be roughly the next 7180 measurements (2592000/361). My record is 
% 30 days long.
date0=datenum(1970,1,1); % give reference date (first date)
kAugustStart = floor((7/12)*82237); 
kAugustEnd = kAugustStart + 7180; 
timeAugust = time(kAugustStart:kAugustEnd, 1);
timeAugust = double(timeAugust)/86400+date0; % in units of days (conversion: seconds/day)
%timeAugustHours = double(timeAugust)/3600 + date0;
pressureAugust = pressure(kAugustStart:kAugustEnd, 1);
% 7181 points in record at 361 second intervals 
%total duration is 7181*361/(24*3600) days
frequency = (0:7180)/(7181*361)*(24*3600); 

% Fourier transforming my data

figure
%plot(pressureAugust);
f = fft(pressureAugust); % Fourier transform of my data
%hold on;
%plot(f) % plots real vs. imaginary parts of data
%plot(abs(f))
semilogy(frequency,abs(real(f)),'-b')
hold on
semilogy(frequency,abs(imag(f)),'-m');
xlabel('cycles per day (cpd)');
xlim([0 125]);
%plot(1/timeAugust,real(f),'-b',1/timeAugust,imag(f),'-m');
%plot(abs(f));
%x_new = ifft(f); % inverse of the Fourier transform

% peaks I can find are at 1 cpd and 2 cpd
% what frequencies correspond to these peaks? - 1 cpd, 2 cpd
% Yes, they're what I'd expect based on the known tidal frequencies, which
% are ~1 cpd for the O1 and K1 tides, and ~2 cpd for the principal lunar
% tide.

%% mean pressure and amplitudes of major peaks

% using Fourier coefficients to find mean pressure, amplitudes of peaks









