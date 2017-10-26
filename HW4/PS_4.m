% SIOC 221A: Analysis of Physical Oceanographic Data Homework #4
% created: Monday, 23 October 2017
% Author: Annie Adelson
% Collaborators: discussed and debugged with


%% data import
%read in automatic data
%initialize vectors
time = [];
pressure = [];
pressure_f = [];

%%read in data

filename = strcat('/Users/AnneElizabeth/Documents/MATLAB/SIO_221A/PS_3/scripps_pier-2015.nc');
time = vertcat(time, ncread(filename,'time'));
pressure = double(vertcat(pressure, ncread(filename,'pressure')));
pressure_f = vertcat(pressure_f, ncread(filename,'pressure_flagPrimary'));



% Define reference date
% Time is measured in seconds since January 1, 1970


%% Filter data

 k = find(pressure_f == 1);
 p = pressure(k);
 t = double(time(k));
 date0=datenum([1970 1 1 0 0 0]);
 date0 = double(date0);
 tplot = (t/24/3600)+date0;
 
 %% time differences
t = double(t);
t_diff = diff(t);
t_diff2 = t_diff/60;
big_diff = max(t_diff2);

%% consider a period with equal increments

cutoff = find(t_diff2>t_diff2(1),1);
p_cut = p(1:cutoff-1);
t_cut = t(1:cutoff-1);
tplot_cut = t_cut/24/3600+date0;


figure 
hold on

plot(tplot_cut, p_cut, 'LineWidth',1)

ax = gca;
xlabel('\fontsize{12} 2015')
ylabel('\fontsize{12}Pressure (dbar)')
title('\fontsize{14}2015 Scripps Pier Sea Water Pressure')
datetick('x','mmm dd','keepticks')
 

%% Fourier transform your data

% frequency domain
dt= t_diff(1);
N=length(p_cut);
T=(dt*N)/(24*3600); % TOTAL time in days
df=1/T;
fn=1/(2*dt);
freq =(0:N/2-1)/(N*dt)*(24*3600);

fft_data = fft(p_cut);
fft_data = fft_data/N;
fft_ddata = fft(diff(p_cut));
fft_ddata = fft_ddata/N;
cutoff = floor(length(p_cut)/2);

figure
subplot(3,1,1)
semilogy(freq(1:cutoff),abs(real(fft_data(1:cutoff))));
xlabel('\fontsize{12} frequency (cycles per day)')
ylabel('\fontsize{12} Pressure (dbar)')
title('\fontsize{14} Real component of the Fourier Transformed Scripps Pier pressure data')
hold on
subplot(3,1,2)
semilogy(freq(1:cutoff),abs(imag(fft_data(1:cutoff))));
xlabel('\fontsize{12} frequency (cycles per day)')
ylabel('\fontsize{12} Pressure (dbar)')
title('\fontsize{14} Imaginary component of the Fourier Transformed Scripps Pier pressure data')
subplot(3,1,3)
semilogy(freq(1:cutoff),abs(fft_data(1:cutoff)));

xlabel('\fontsize{12} frequency (cycles per day)')
ylabel('\fontsize{12} Pressure (dbar)')
title('\fontsize{14} Absolute value of the Fourier Transformed Scripps Pier pressure data')

o1 = 1/25.82*24;
k1 = 1/23.93*24;
m2 = 1/12.42*24;

liney = 1e-12:1e6:1e10;



figure
semilogy(freq(1:cutoff),abs(fft_data(1:cutoff)));
hold on


pks = abs(fft_data(1:cutoff));
[Ma,I] = max(pks); %largest amplitude
[Mb,Ib] = max(pks(pks~=Ma));
[Mc,Ic]= max(pks(pks~=Ma & pks ~= Mb));
[Md,Id]= max(pks(pks~=Ma & pks ~= Mb & pks ~= Mc));

plot(k1*ones(size(liney)),liney,'-k'); 
plot(o1*ones(size(liney)),liney,'-.k');
plot(m2*ones(size(liney)),liney,':k');

plot(freq(find(abs(fft_data) == Ma,1)), Ma, 'or');
plot(freq(find(abs(fft_data) == Mb,1)), Mb, 'or');
plot(freq(find(abs(fft_data) == Mc,1)), Mc, 'or');
plot(freq(find(abs(fft_data) == Md,1)), Md, 'or');



%[pks, locs] = findpeaks(abs(fft_data(1:cutoff)));
%[M,I] = max(pks); %largest amplitude
%[Mb,Ib] = max(pks(pks~=max(pks)));
%[Mc,Ic]= max(pks(pks~=max(pks) & pks ~= Mb));

xlabel('\fontsize{12} frequency (cycles per day)')
ylabel('\fontsize{12} Pressure (dbar)')
legend('Fourier transformed data', 'K1 frequency', 'O1 frequency', 'M2 frequency','peaks')
xlim([0 3])
ylim([1e-4 1e1])

title('\fontsize{14}2015 Scripps Pier Sea Water Pressure')


%% normalized spectral density 
figure
loglog(freq(1:cutoff),abs(fft_data(1:cutoff)).^2);
hold on


xlabel('\fontsize{12} frequency (cycles per day)')
ylabel('\fontsize{12} Spectral density')
legend('Fourier transformed data', 'peaks')

title('\fontsize{14}2015 Scripps Pier Sea Water Pressure')


%% compare to least squares fit
%o1_amp_winter = 0.168915728632706
%k1_amp_winter = 0.371008365460424
%m2_amp_winter = 0.527422224875983

o1_peak = 2*Md;
k1_peak = 2*Mc;
m2_peak = 2*Mb;

%% ddata

figure
semilogy(freq(1:cutoff),abs(fft_ddata(1:cutoff)));


