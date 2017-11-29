% file SIOC 221A HW 8
% 
% author Julia Dohner, with help from Annie Adelson
%
% due date November 30, 2017

clear all; close all;

% TODO: whether to put if statement before or after alias period
% calculation


%% aliasing


f_M2=24/12.42; % cycles per day
f_S1 = 24/24;

% fast sampling
f_sampling_fast=1/(0.99349); % cycles per day
f_Nyquist_fast=f_sampling_fast/2;

% M2 aliasing in fast sampling
M_fast_M2=floor(f_M2/f_Nyquist_fast); % compute the integer ratio of the two frequencies.
alias_fast_M2 = f_M2 - floor(f_M2/f_Nyquist_fast)*f_Nyquist_fast;
% Note: if M is odd then reset - What is this?? Why not before above line?
if(rem(M_fast_M2,2)~=0) 
    alias_fast_M2=f_Nyquist_fast-alias_fast_M2; 
end
alias_period_fast_M2 = 1/alias_fast_M2

% S1 aliasing in fast sampling
M_fast_S1=floor(f_S1/f_Nyquist_fast); % compute the integer ratio of the two frequencies.
alias_fast_S1 = f_S1 - floor(f_S1/f_Nyquist_fast)*f_Nyquist_fast;
% Note: if M is odd then reset
if(rem(M_fast_S1,2)~=0) 
    alias_fast_S1=f_Nyquist_fast - alias_fast_S1; 
end
alias_period_fast_S1 = 1/alias_fast_S1


% Fast sampling, M2: M = 3, alias_period = 2.3667 days
% Fast sampling, S1: M = 1, alias_period = 2.0132


% science sampling
f_sampling_sci=1/(20.86455); % cycles per day
f_Nyquist_sci=f_sampling_sci/2;

% M2 aliasing in science sampling
M_sci_M2=floor(f_M2/f_Nyquist_sci); % compute the integer ratio of the two frequencies.
alias_sci_M2 = f_M2 - floor(f_M2/f_Nyquist_sci)*f_Nyquist_sci;
% Note: if M is odd then reset
if(rem(M_sci_M2,2)~=0) 
    alias_sci_M2=f_Nyquist_sci - alias_sci_M2; 
end
alias_period_sci_M2 = 1/alias_sci_M2


% S1 aliasing in science sampling
M_sci_S1=floor(f_S1/f_Nyquist_sci); % compute the integer ratio of the two frequencies.
alias_sci_S1 = f_S1 - floor(f_S1/f_Nyquist_sci)*f_Nyquist_sci;
% Note: if M is odd then reset
if(rem(M_sci_S1,2)~=0) 
    alias_sci_S1=f_Nyquist_sci - alias_sci_S1; 
end
alias_period_sci_S1 = 1/alias_sci_S1

%This calculation shows that M is 80 and the tidal energy aliases into a
% 65.6178 day period (alias_period)


%% Frequency-wavenumber spectra

time = [];
lat = [];
lon = [];
sst = [];

% time units = 'days since 1800-1-1 00:00:00'
time = [time; ncread(strcat('sst.mnmean.nc'),'time')];
% units are degrees north
lat = [lat; ncread(strcat('sst.mnmean.nc'),'lat')];
% units are degrees east
lon = [lon; ncread(strcat('sst.mnmean.nc'),'lon')];
% sst dimensions are lon,lat,time (360x180x415)
sst = [sst; ncread(strcat('sst.mnmean.nc'),'sst')];

date0=datenum(1800,1,1); % give reference date (first date)
time2 = double(time)+date0; 

% Compute frequency/wavenumber spectra for two slabs of data in the Pacific.
% Along the equator: -0.5?S, from 219.5 to 269.5?E.
% In the subtropics: -29.5?S, from 189.5 to 249.5?E.

% For both slabs, plot a time/longitude plot (known as a Hovmoller diagram)
% and compute frequency/wavenumber spectra using at least some amount of 
% averaging/segmenting.

% slab a: Along the equator: -0.5?S, from 219.5 to 269.5?E.
lat_a_ind = find(lat == -0.5);
lon_a_ind = find(lon == 219.5);
lat_a = lat(lat_a_ind);
x = 219.5:269.5; % calculating how far along to go in lon vector
last_lon_a = lon_a_ind + length(x) -1;
lon_a = lon(lon_a_ind:last_lon_a);
% sst dimensions are lon,lat,time (360x180x415)
time_a = time2(:);
sst_a = sst(lon_a_ind:last_lon_a,lat_a_ind,:);



sst_a2 = squeeze(sst_a);
sst_a3 = sst_a2;
% time on y axis, location in x, sst in color
figure('name','Along the equator: -0.5?S, from 219.5 to 269.5?E.');
pcolor(time_a,lon_a,sst_a3);
shading interp;
colorbar;
set(gca,'FontSize',16);
title('Along the equator: -0.5?S, from 219.5 to 269.5?E');
xlabel('time');
ylabel('Longitude');
datetick('x','mm/dd/yy')
axis tight

% slab b: In the subtropics: -29.5?S, from 189.5 to 249.5?E.
lat_b_ind = find(lat == -29.5);
lon_b_start = find(lon == 189.5);
lat_b = lat(lat_b_ind);
y = 189.5:249.5; % calculating how far along to go in lon vector
last_lon_b = lon_b_start + length(y) -1;
lon_b = lon(lon_b_start:last_lon_b);
% sst dimensions are lon,lat,time (360x180x415)
time_b = time2(:);
sst_b = sst(lon_b_start:last_lon_b,lat_b_ind,:);

sst_b2 = squeeze(sst_b);
sst_b3 = sst_b2;
% time on y axis, location in x, sst in color
figure('name','In the subtropics: -29.5?S, from 189.5 to 249.5?E.');

pcolor(time_b,lon_b,sst_b3);
shading interp;
colorbar;
set(gca,'FontSize',16);
title('In the subtropics: -29.5?S, from 189.5 to 249.5?E.');
xlabel('time');
ylabel('Longitude');
datetick('x','mm/dd/yy')
axis tight

% question: do we want time on x or y axis?

%% 
% compute frequency/wavenumber spectra using at least some amount of 
% averaging/segmenting


% slab a: Along the equator: -0.5?S, from 219.5 to 269.5?E.
% pcolor(time_a,lon_a,sst_a3);
% looking at one latitude, so just 3D (lon, time, sst)
% looking at one property (sst) over two dimensions
% axes are longitude (wavenumber) and time (frequency)
% need to convert longitude and time to wavenumber and frequency
% sst_a3 is currently 
% segmenting in time - this is the 415 dimension

N_data = floor(length(time2)/7)*7; % use some number of datapoints divisible by 8
N = N_data/7; % length of each chunk of data (aka segment length)
M = N_data/N; % number of segments splitting data into

% break data into segments
sst_a4 = sst_a3(1:end-1,1:N_data);
sst_a5 = reshape(sst_a4, 50, N, M);
ft_a_time = fft(sst_a5);
%ft_a_space = fft(sst_a5,
% question: how do we do the fft in space? which dimension?
% compute squared amplitude for half of fft
a_time_amp = (abs(ft_a_time(1:floor(length(ft_a_time)/2)+1,:)).^2); % +1 for even N
% Matthew doesn't specify fft

%force them to be even - makes everything easier
% if rem(length(iz),2)==1
%     iz=iz(1:end-1);
% end


