% file SIOC 221A HW 9
% 
% author Julia Dohner
%
% due date December 13, 2017

% analysis of weekly flask-collected co2 data at Mauna Loa and La Jolla
% stations

clear all; close all;



%% load CO2 data

dataMLO = fopen('Early_LaJolla_CO2_halfhourly_1958-1962.csv');

% read in header info
headerInfo = textscan(dataMLO, '%s', 37, 'delimiter', '\n');

% read in the data
valsMLO = textscan(dataMLO, '%f %f %f %f %f %f', ...
    'delimiter',',,,,,');


fclose(dataMLO);

% format of .txt files is year, co2 value
MLOyear = valsMLO{1};
MLOmonth = valsMLO{2};
MLOday = valsMLO{3};
MLOhour = valsMLO{4};
MLOminute = valsMLO{5};
MLOco2 = valsMLO{6};

% cut off data to start at 1990 because dataset too big
startIndex = find(MLOyear == 2004 & MLOmonth == 12 & MLOday == 29,1);
MLOco2_2 = MLOco2(startIndex:end);
MLOyear_2 = MLOyear(startIndex:end);
MLOmonth_2 = MLOmonth(startIndex:end);
MLOday_2 = MLOday(startIndex:end);
MLOhour_2 = MLOhour(startIndex:end);
MLOminute_2 = MLOminute(startIndex:end);
MLOseconds = zeros(length(MLOyear_2),1); % placeholder for datenum since no seconds data

% remove nan's
addpath('/Users/juliadohner/Documents/MATLAB/SIOC_221A/HW9/Inpaint_nans/Inpaint_nans');
MLOco2_2 = inpaint_nans(MLOco2_2);

%converting YMD into Matlab's time interpretation 
dateData = datenum(MLOyear_2, MLOmonth_2, MLOday_2, MLOhour_2, MLOminute_2, MLOseconds);


%% inspect spacing of data

% inspect time spacing between measurements
MLO_t_diff = diff(MLOminute);

% figure
%plot(1:length(minute)-1,MLO_t_diff);

minDiffMLO = min(MLO_t_diff);
maxDiffMLO = max(MLO_t_diff);


% plot timeseries
figure('name','Atmospheric CO2 Timeseries');
plot(dateData,MLOco2_2,'-b');
xlabel('\fontsize{14}year')
ylabel('\fontsize{14}ppm')
title('\fontsize{16}Atmospheric CO2 Records')
legend('\fontsize{12}La Jolla Station','\fontsize{12}Mauna Loa', '\fontsize{12}South Pole','Location','northwest');


%% compute spectra (from in-class coherence example)

% three overlapping segments:
N = length(LJOco2_2);
Nseg = 4; % number of segments splitting data into
segment_length = N/Nseg; % length of each chunk of data (aka segment length)
M = segment_length/2;

LJO_use = [reshape(LJOco2_2,segment_length,Nseg)];
MLO_use=[reshape(MLOco2_2,segment_length,Nseg)];
SPO_use=[reshape(SPOco2_2,segment_length,Nseg)];

LJO_ft=fft(detrend(LJO_use).*(hann(segment_length)*ones(1,Nseg)));
MLO_ft=fft(detrend(MLO_use).*(hann(segment_length)*ones(1,Nseg)));
SPO_ft=fft(detrend(SPO_use).*(hann(segment_length)*ones(1,Nseg)));


% % using three segments with 50% overlap (concatenating in lines below)
% LJO_use=[reshape(LJOco2_2,segment_length,Nseg) reshape(LJOco2_2(M+1:end-M),segment_length,Nseg-1)];
% MLO_use=[reshape(MLOco2_2,segment_length,Nseg) reshape(MLOco2_2(M+1:end-M),segment_length,Nseg-1)];
% SPO_use=[reshape(SPOco2_2,segment_length,Nseg) reshape(SPOco2_2(M+1:end-M),segment_length,Nseg-1)];
% 
% LJO_ft=fft(detrend(LJO_use).*(hann(segment_length)*ones(1,3)));
% MLO_ft=fft(detrend(MLO_use).*(hann(segment_length)*ones(1,3)));
% SPO_ft=fft(detrend(SPO_use).*(hann(segment_length)*ones(1,3)));

% question: mean is still in here. Do I want it?
LJO_spec=sum(abs(LJO_ft(1:M+1,:)).^2,2)/N; % sum over all spectra (2nd dim)
LJO_spec(2:end)=LJO_spec(2:end)*2; % multiply by 2 to make up for lost energy

MLO_spec=sum(abs(MLO_ft(1:M+1,:)).^2,2)/N;
MLO_spec(2:end)=MLO_spec(2:end)*2;

SPO_spec=sum(abs(SPO_ft(1:M+1,:)).^2,2)/N;
SPO_spec(2:end)=SPO_spec(2:end)*2;

%% uncertainty estimate on spectra

% question: is this nu calculation right??
nu = 1.9*2*Nseg; % DOF = 2*number of segments*1.9 (1.9 for the Hanning)
err_high = nu/chi2inv(0.05/2,nu);
err_low = nu/chi2inv(1-0.05/2,nu);
ratio_chi2 = err_high/err_low;

%% plot spectra

% Question: I added the +1 to the M here to make it match the dims of the specs
% (145x1) because means are still in there
frequency=(1:M+1)/(segment_length/12);
frequency = frequency';

figure('name','Power Spectra of CO2 Records');
% TODO: fix the locations of the uncertainty estimates
loglog(frequency,LJO_spec, '-r', [.2 .2],[err_low err_high]*LJO_spec(50), ...
    frequency, MLO_spec, '-b',[.1 .1],[err_low err_high]*MLO_spec(50), ...
    frequency, SPO_spec, '-g',[.3 .3],[err_low err_high]*SPO_spec(50))
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}ppm^2/cpy')
title('\fontsize{16}Power Spectra of CO2 Records')
legend('\fontsize{12}La Jolla Station','\fontsize{12}Mauna Loa Station', '\fontsize{12}South Pole','Location','northeast');

%% compute coherence

% compute cross covariance of LJO/MLO
ccLM=sum(LJO_ft(1:M+1,:).*conj(MLO_ft(1:M+1,:)),2)/N;
ccLM(2:end)=ccLM(2:end)*2;
% compute coherence
C_LM=abs(ccLM)./sqrt(LJO_spec.*MLO_spec);
phase_LM = atan2(-imag(ccLM),real(ccLM));
deltaPhase_LM = sqrt((1-C_LM.^2)./(abs(C_LM).^2*2*Nseg)); 

% compute cross covariance of MLO/SPO
ccMS=sum(MLO_ft(1:M+1,:).*conj(SPO_ft(1:M+1,:)),2)/N;
ccMS(2:end)=ccMS(2:end)*2;
% compute coherence
C_MS=abs(ccMS)./sqrt(MLO_spec.*SPO_spec);
phase_MS = atan2(-imag(ccMS),real(ccMS));
deltaPhase_MS = sqrt((1-C_MS.^2)./(abs(C_MS).^2*2*Nseg));

%% uncertainty estimate 

alpha = 0.05; % 95% confidence level
gamma_threshold= sqrt(1-alpha^(1/(Nseg-1)));


%% plot the coherence

figure('name','Coherence Plots of CO2 Records');
subplot(2,1,1)
plot(frequency, C_LM,[frequency(1) frequency(end)],[gamma_threshold gamma_threshold]);
axis tight;
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}coherence')
title('\fontsize{16}Coherence of CO2 Records (LJO and MLO)')
legend('\fontsize{12}La Jolla vs. Mauna Loa','\fontsize{12}Uncertainty Threshold','Location','southeast');

subplot(2,1,2)
plot(frequency, C_MS,[frequency(1) frequency(end)],[gamma_threshold gamma_threshold]); 
axis tight
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}coherence')
title('\fontsize{16}Coherence of CO2 Records (MLO and SPO)')
legend('\fontsize{12}Mauna Loa vs. South Pole','\fontsize{12}Uncertainty Threshold','Location','southeast');

% Question: is this threshold right? seems really high... but I suppose it's that 5%
% of the data will be above the threshold just due to random chance, so it should be
% high?
% "set a threshold for evaluating whether a calculated coherence exceeds 
% what we might expect from random white noise"

%% plot the phase

% Question: I'm having a hard time understanding this plot (read Lec 15
% Notes)
% Question: what are the units on the y axis? (again, read lec 15 notes)
figure('name','Phase Plot for Mauna Loa and South Pole Coherence');
semilogx(frequency,[phase_MS phase_MS+deltaPhase_MS phase_MS-deltaPhase_MS]);
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}phase units')
title('\fontsize{16}Phase Plot for Mauna Loa and South Pole Coherence')
legend('\fontsize{12}phase', 'phase+delta phase','phase-delta phase','Location','northeast');
axis tight