% file SIOC 221A HW 8
% 
% author Julia Dohner, with help from Annie Adelson
%
% due date November 30, 2017

clear all; close all;


%% aliasing


f_M2=24/12.42; % cycles per day
f_S1 = 24/24;

% fast sampling
f_sampling_fast=1/(0.99349); % cycles per day
f_Nyquist_fast=f_sampling_fast/2;

% M2 aliasing in fast sampling
M_fast_M2=floor(f_M2/f_Nyquist_fast); % compute the integer ratio of the two frequencies.
alias_fast_M2 = f_M2 - floor(f_M2/f_Nyquist_fast)*f_Nyquist_fast;
% Note: if M is odd then reset 
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

% slab a: Along the equator: -0.5 S, from 219.5 to 269.5 E.
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
figure('name','Along the equator: -0.5\circS, from 219.5 to 269.5\circE');
pcolor(time_a,lon_a,sst_a3);
shading interp;
h = colorbar; ylabel(h,'\fontsize{14}\circCelsius');
set(gca,'FontSize',16);
title('\fontsize{16}Along the equator: -0.5\circS, from 219.5 to 269.5\circE');
xlabel('\fontsize{14}time');
ylabel('\fontsize{14}degrees longitude');
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
figure('name','In the subtropics: -29.5\circS, from 189.5 to 249.5\circE');
pcolor(time_b,lon_b,sst_b3);
shading interp;
h = colorbar; ylabel(h,'\fontsize{14}\circCelsius');
set(gca,'FontSize',16);
title('\fontsize{16}In the subtropics: -29.5\circS, from 189.5 to 249.5\circE.');
xlabel('\fontsize{14}time');
ylabel('\fontsize{14}degrees longitude');
datetick('x','mm/dd/yy')
axis tight


%% 
% compute frequency/wavenumber spectra using at least some amount of 
% averaging/segmenting


% slab a: Along the equator: -0.5 S, from 219.5 to 269.5 E.

N_data = floor(length(time2)/6)*6; % use some number of datapoints divisible by 6
N = N_data/6; % length of each chunk of data (aka segment length)
M = N_data/N; % number of segments splitting data into
p = 50; % number of rows (in this case number of longitude values)

% break data into segments
sst_a4 = sst_a3(1:end-1,1:N_data);
sst_a5 = reshape(sst_a4, p, N, M);

% compute fft in time (inner fft) and space (outer fft)
st_a = fftshift(fft2(sst_a5))./p./N; % fft in 2-D, normalizing

% time and space increments
t_diff = diff(time);
t_diff_mean = mean(t_diff);
dt = 1/t_diff_mean; % time interval
dz = 1; % 1 degree longitude

% fundamental frequency and wavenumber
%df=1./N./dt; % time 
df = 1/N/dt*dt*12; % in cycles per year
dk=1./p./dz; % space

% average amplitudes for all realizations 
% amplitudes: (*conj is same as abs of value^2)
amp_a=st_a.*conj(st_a)./df./dk;

% creating frequency and wavenumber vectors for plotting purposes
f=[-fliplr(1:(N/2)) 0 (1:(N/2-1))].*df; % frequency
k=[-fliplr(1:(p/2)) 0 (1:(p/2-1))].'.*dk;

spec_a_avg = mean(amp_a,3); % average amplitudes
figure('name','Frequency-wavenumber spectrum along the equator (-0.5 S, from 219.5 to 269.5 E)');
imagesc(f,k,log(spec_a_avg)); % y axis is longitude, x axis is time
h = colorbar;
title('\fontsize{16}Frequency-wavenumber spectrum along the equator (-0.5\circS, from 219.5 to 269.5\circE)');
xlabel('\fontsize{14}frequency (cycles per year)'); % data taken monthly
ylabel('\fontsize{14}wavenumber (cycles per degree longitude)');

%% slab b: In the subtropics: -29.5 S, from 189.5 to 249.5 E.

N_data = floor(length(time2)/6)*6; % use some number of datapoints divisible by 6
N = N_data/6; % length of each chunk of data (aka segment length)
M = N_data/N; % number of segments splitting data into
p_b = 60; % number of rows (in this case number of longitude values)

% break data into segments
sst_b4 = sst_b3(1:end-1,1:N_data);
sst_b5 = reshape(sst_b4, p_b, N, M);

% compute fft in time (inner fft) and space (outer fft)
st_b = fftshift(fft2(sst_b5))./p_b./N; % fft in 2-D, normalizing

% time and space increments
t_diff = diff(time);
t_diff_mean = mean(t_diff);
dt = 1/t_diff_mean; % time interval
dz = 1; % 1 degree longitude

% fundamental frequency and wavenumber
%df=1./N./dt; % time 
df = 1/N/dt*dt*12; % in cycles per year
dk=1./p./dz; % space

% average amplitudes for all realizations 
% amplitudes: (*conj is same as abs of value^2)
amp_b=st_b.*conj(st_b)./df./dk;

% creating frequency and wavenumber vectors for plotting purposes
f=[-fliplr(1:(N/2)) 0 (1:(N/2-1))].*df; % frequency
k=[-fliplr(1:(p_b/2)) 0 (1:(p_b/2-1))].'.*dk;

spec_b_avg = mean(amp_b,3); % average amplitudes
figure('name','Frequency-wavenumber spectrum in the subtropics (-29.5 S, from 189.5 to 249.5 E)');
imagesc(f,k,log(spec_b_avg)); % y axis is longitude, x axis is time
h = colorbar;
title('\fontsize{16}Frequency-wavenumber spectrum in the subtropics: (-29.5\circS, from 189.5 to 249.5\circE)');
xlabel('\fontsize{14}frequency (cycles per year)'); % data taken monthly
ylabel('\fontsize{14}wavenumber (cycles per degree longitude)');