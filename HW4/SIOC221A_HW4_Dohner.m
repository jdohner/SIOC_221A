% file SIOC 221A HW 4
% 
% author Julia Dohner
%
% due date October 26, 2017

clear all; close all;

% TODO: 
% half the pressure vector from the beginning rather than only plotting first half
% 1. get third peak
% 2. Having trouble with getting amplitudes via indexing - why do my peaks b and d have different amplitudes but the
% same index?
% 3. other way to see if spectral peaks align with LSF?



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
% total duration is 7181*361/(24*3600) days


%% Fourier transforming my data


f = fft(pressureAugust); % Fourier transform of my data
data = f(1:3590);
N = length(pressureAugust);
frequency = (0:7180)/(7181*361)*(24*3600); 
frequency = frequency(1:3590)'; % only using first half of vector
data = data/N; % Normalizing the data


figure
semilogy(frequency,abs(real(data)),'-b')
hold on
semilogy(frequency,abs(imag(data)),'-m');
title('Real and imaginary parts of Fourier transform of Scripps Pier pressure in August 2015');
xlabel('cycles per day (cpd)');
ylabel('dbar');
%xlim([0 125]);

figure
semilogy(frequency,abs(data),'-b')
title('Fourier transform of Scripps Pier pressure in August 2015');
xlabel('cycles per day (cpd)');
ylabel('dbar');
xlim([0 125]);

% question: how to resolve the third peak??

% peaks I can find are at 1 cpd and 2 cpd
% what frequencies correspond to these peaks? - 1 cpd, 2 cpd
% Yes, they're what I'd expect based on the known tidal frequencies, which
% are ~1 cpd for the O1 and K1 tides, and ~2 cpd for the principal lunar
% tide.

% sampling frequency = 1/timestep



%% mean pressure and amplitudes of major peaks
% using Fourier coefficients to find mean pressure, amplitudes of peaks

% mean pressure is the first value in the Fourier coefficient vector
meanPressure = data(1);

absData = abs(data);

[amp ind] = sort(absData,'descend');
amplitudes(:,2) = amp(1:4);
amplitudes(:,1) = ind(1:4); % stores the indices (col 1) and amplitudes (col 2) of major peaks

% multiple amplitudes by 2 to account for both positive and negative
% frequencies
amplitude_a = 2*amplitudes(1,2); % this is the mean
amplitude_b = 2*amplitudes(2,2); % largest non-mean peak
amplitude_c = 2*amplitudes(3,2); % second-largest non-mean peak
amplitude_d = 2*amplitudes(4,2); % third-largest non-mean peak

% TODO: my aplitudes seem to be too small....


%% Spectral peaks vs. least-squares fit?

% They align in terms of their frequency and amplitude. Otherwise I'm not
% sure what I could do to further compare them.


%% spectral energy vs. frequency

spectralE = abs(data).^2;

figure
semilogy(frequency,spectralE,'-b')
title('Spectral Energy of Scripps Pier pressure in August 2015');
xlabel('cycles per day (cpd)');
ylabel('dbar^2 / cpd');
xlim([0 125]);

% Question: doesn't look much different from my FT plot

% the spectrum is red (peaks at low frequencies)

% differentiating time series in time before FT
% time needs to be a positive integer - how to do this? convert back to
% seconds?

timeAugust_diff = time(kAugustStart:kAugustEnd, 1);
timeAugust_diff = double(timeAugust)+date0; % in units of seconds % ROOM FOR BUG HERE

% getting error
% Error using diff
% Difference order N must be a positive integer scalar.
diff_august = diff(pressureAugust);
f_diff = fft(diff_august); % Fourier transform of my data
data_diff = f_diff(1:3590);
data_diff = data_diff/N; % Normalizing the data

spectralE_diff = abs(data_diff).^2;

% TODO: unsure of axis labels for the first derivative

figure
semilogy(frequency,spectralE_diff,'-b')
title('Spectral Energy of time derivative of Scripps Pier pressure in August 2015');
xlabel('cycles per day (cpd)');
ylabel('(dbar/time)^2 / cpd');
xlim([0 125]);
