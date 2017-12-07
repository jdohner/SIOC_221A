% file SIOC 221A HW 9
% 
% author Julia Dohner
%
% due date December 13, 2017

% analysis of weekly flask-collected co2 data at Mauna Loa and La Jolla
% stations

clear all; close all;

dataMLO = fopen('monthly_flask_co2_mlo_JLD.txt');
dataLJO = fopen('monthly_flask_co2_ljo_JLD.txt');
dataSPO = fopen('monthly_flask_co2_spo_JLD.txt');

valsMLO = textscan(dataMLO, '%f %f', ...
    'delimiter','\t');
valsLJO = textscan(dataLJO, '%f %f', ...
    'delimiter','\t');
valsSPO = textscan(dataSPO, '%f %f', ...
    'delimiter','\t');

fclose(dataMLO);
fclose(dataLJO);
fclose(dataSPO);

% format of .txt files is year, flag, co2 value
LJOyear = valsLJO{1};
LJOco2 = valsLJO{2};

MLOyear = valsMLO{1};
MLOco2 = valsMLO{2};

SPOyear = valsSPO{1};
SPOco2 = valsSPO{2};

% remove flagged data
% defunct now that I've chosen data with no NaN's
for i = 1:length(MLOco2)
    if MLOco2(i) == -99.99
        MLOco2(i) = nan;
    end
end
for i = 1:length(LJOco2)
    if LJOco2(i) == -99.99
        LJOco2(i) = nan;
    end
end
for i = 1:length(SPOco2)
    if SPOco2(i) == -99.99
        SPOco2(i) = nan;
    end
end

% remove nan's
% TODO: if time, go back and change this to a linear interpolation
addpath('/Users/juliadohner/Documents/MATLAB/SIOC_221A/HW9/Inpaint_nans/Inpaint_nans');
MLOco2 = inpaint_nans(MLOco2);
LJOco2 = inpaint_nans(LJOco2);
SPOco2 = inpaint_nans(SPOco2);

% plot timeseries
figure
plot(LJOyear,LJOco2, '-r', MLOyear,MLOco2,'-b',SPOyear,SPOco2,'-g');
xlabel('\fontsize{14}year')
ylabel('\fontsize{14}ppm')

% inspect time spacing between measurements
MLO_t_diff = diff(MLOyear);
LJO_t_diff = diff(LJOyear);
SPO_t_diff = diff(SPOyear);

% figure
% plot(1:length(MLOyear)-1,MLO_t_diff);
% figure
% plot(1:length(LJOyear)-1,LJO_t_diff);
% figure
% plot(1:length(SPOyear)-1,SPO_t_diff);

minDiffLJO = min(LJO_t_diff);
maxDiffLJO = max(LJO_t_diff);
minDiffMLO = min(MLO_t_diff);
maxDiffMLO = max(MLO_t_diff);
minDiffSPO = min(SPO_t_diff);
maxDiffSPO = max(SPO_t_diff);

% too big differences in time increments of daily LJO flask sampling, will try 
% monthly data instead

% update: time differences range between 0.0767 and 0.0850 years for both 
% MLO and LJO records (difference of 0.0083) or ~3 days, deemed small 
% enough difference in time increments
% read: counts in my book as even spacing

startYear = LJOyear(1,1); % 1.969041100000000e+03 (latest start)
endYear = SPOyear(length(SPOyear),1); % earliest end
startIndex_LJO = find(LJOyear == startYear);
startIndex_MLO = find(MLOyear == startYear);
startIndex_SPO = find(SPOyear == startYear);
endIndex_LJO = find(LJOyear == endYear);
endIndex_MLO = find(MLOyear == endYear);
endIndex_SPO = find(SPOyear == endYear);

% create new vectors
LJOco2_2 = LJOco2(startIndex_LJO:endIndex_LJO);
LJOyear_2 = LJOyear(startIndex_LJO:endIndex_LJO);
MLOco2_2 = MLOco2(startIndex_MLO:endIndex_MLO);
MLOyear_2 = MLOyear(startIndex_MLO:endIndex_MLO);
SPOco2_2 = SPOco2(startIndex_SPO:endIndex_SPO);
SPOyear_2 = SPOyear(startIndex_SPO:endIndex_SPO);


%% compute spectra (from in-class coherence example)

N = length(LJOco2_2);
segment_length = N/2; % length of each chunk of data (aka segment length)
M = segment_length/2;
Nseg = N/segment_length; % number of segments splitting data into

% using three segments with 50% overlap (concatenating in lines below)
LJO_use=[reshape(LJOco2_2,segment_length,Nseg) reshape(LJOco2_2(M+1:end-M),segment_length,Nseg-1)];
MLO_use=[reshape(MLOco2_2,segment_length,Nseg) reshape(MLOco2_2(M+1:end-M),segment_length,Nseg-1)];
SPO_use=[reshape(SPOco2_2,segment_length,Nseg) reshape(SPOco2_2(M+1:end-M),segment_length,Nseg-1)];

LJO_ft=fft(detrend(LJO_use).*(hann(segment_length)*ones(1,3)));
MLO_ft=fft(detrend(MLO_use).*(hann(segment_length)*ones(1,3)));
SPO_ft=fft(detrend(SPO_use).*(hann(segment_length)*ones(1,3)));

% note: mean is still in here. Do I want it?
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

% I added the +1 to the M here to make it match the dims of the specs
% (145x1) because means are still in there
frequency=(1:M+1)/(segment_length/12);
frequency = frequency';

% see harmonics on the annual cycle

figure('name','spectra with 3 segments, 50% overlap');
loglog(frequency,LJO_spec, '-r', [.2 .2],[err_low err_high]*LJO_spec(100), ...
    frequency, MLO_spec, '-b',[.1 .1],[err_low err_high]*MLO_spec(100), ...
    frequency, SPO_spec, '-g',[.3 .3],[err_low err_high]*SPO_spec(100))
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}ppm^2/cpy')
title('\fontsize{16}Spectrum of Mauna Loa and La Jolla CO2 records')
legend('\fontsize{12}La Jolla Station','\fontsize{12}Mauna Loa Station', '\fontsize{12}South Pole Observatory');

%% compute coherence

% compute cross covariance of LJO/MLO
ccLM=sum(LJO_ft(1:M+1,:).*conj(MLO_ft(1:M+1,:)),2)/N;
ccLM(2:end)=ccLM(2:end)*2;
% compute coherence
C_LM=abs(ccLM)./sqrt(LJO_spec.*MLO_spec);
phase_LM = atan2(-imag(ccLM),real(ccLM));
deltaPhase_MS = sqrt((1-ccLM.^2)./(abs(ccLM).^2*2*Nseg)); % TODO: having issues with deltaPhase, plotting

% compute cross covariance of MLO/SPO
ccMS=sum(MLO_ft(1:M+1,:).*conj(SPO_ft(1:M+1,:)),2)/N;
ccMS(2:end)=ccMS(2:end)*2;
% compute coherence
C_MS=abs(ccMS)./sqrt(MLO_spec.*SPO_spec);
phase_MS = atan2(-imag(ccMS),real(ccMS));
deltaPhase_MS = sqrt((1-ccMS.^2)./(abs(ccMS).^2*2*Nseg));

%% uncertainty estimate 

alpha = 0.05; % 95% confidence level
gamma_threshold= sqrt(1-alpha^(1/(Nseg-1)));

% TODO: label all stuff on coherence plots
% maybe switch to daily records instead
% need delta on phase?
% interpreting phase: lecture 15 notes
% uncertainty on coherence: also lecture 15 notes



%% plot the coherence

figure
plot(frequency, C_LM);

figure
plot(frequency, C_MS,[0.1 0.6],[gamma_threshold gamma_threshold]); %, frequency,[phase_MS phase_MS+delta_phase_MS phase_MS-delta_phase_MS]);
% having trouble working with delta phase

%% my old code: create 3 segments with 50% overlap


% N = length(LJOco2_2);
% Nseg = 2; % number of segments splitting data into
% segment_length = N/Nseg; % length of each chunk of data (aka segment length)
% M = segment_length/2;
% 
% LJOco2_3 = reshape(LJOco2_2,segment_length,Nseg);
% LJOco2_4 = LJOco2_2(M+1:length(LJOco2_2)-M);
% LJOco2_5 = [LJOco2_3 LJOco2_4];
% LJOco2_5 = detrend(LJOco2_5).*(hann(segment_length)*ones(1,3)); % 3 = number of columns
% LJOco2_ft = fft(LJOco2_5);
% 
% MLOco2_3 = reshape(MLOco2_2,segment_length,Nseg);
% MLOco2_4 = MLOco2_2(M+1:length(MLOco2_2)-M);
% MLOco2_5 = [MLOco2_3 MLOco2_4];
% MLOco2_5 = detrend(MLOco2_5).*(hann(segment_length)*ones(1,3));
% MLOco2_ft = fft(MLOco2_5);
% 
% SPOco2_3 = reshape(SPOco2_2,segment_length,Nseg);
% SPOco2_4 = SPOco2_2(M+1:length(SPOco2_2)-M);
% SPOco2_5 = [SPOco2_3 SPOco2_4];
% SPOco2_5 = detrend(SPOco2_5).*(hann(segment_length)*ones(1,3));
% SPOco2_ft = fft(SPOco2_5);
% 
% % next step: calculate amplitude
% LJOco2_amp = 2*(abs(LJOco2_ft(1:M+1,:)).^2)/segment_length;
% LJOco2_amp = LJOco2_amp(2:length(LJOco2_amp),:); % dump the mean
% LJOco2_mean = mean(LJOco2_amp,2); % take the mean across rows
% 
% MLOco2_amp = 2*(abs(MLOco2_ft(1:M+1,:)).^2)/segment_length;
% MLOco2_amp = MLOco2_amp(2:length(MLOco2_amp),:); % dump the mean
% MLOco2_mean = mean(MLOco2_amp,2); % take the mean across rows
% 
% SPOco2_amp = 2*(abs(SPOco2_ft(1:M+1,:)).^2)/segment_length;
% SPOco2_amp = SPOco2_amp(2:length(SPOco2_amp),:); % dump the mean
% SPOco2_mean = mean(SPOco2_amp,2); % take the mean across rows
% 
% frequency=(1:M)/(segment_length/12);
% frequency = frequency';
% 
% % plot spectra
% % uncertainty estimate
% nu = 2*1; % DOF = 2*number of segments
% % do I need to segment?
% err_high_MLO = nu/chi2inv(0.05/2,nu);
% err_low_MLO = nu/chi2inv(1-0.05/2,nu);
% ratio_chi2_MLO = err_high_MLO/err_low_MLO;
% 
% % see harmonics on the annual cycle
% 
% figure('name','spectra with 3 segments, 50% overlap');
% loglog(frequency,LJOco2_mean, '-r', frequency, MLOco2_mean, '-b',[.1 .1],[err_low_MLO err_high_MLO]*MLOco2_amp(100),frequency, SPOco2_mean, '-g')
% xlabel('\fontsize{14}cycles per year')
% ylabel('\fontsize{14}ppm^2/cpy')
% title('\fontsize{16}Spectrum of Mauna Loa and La Jolla CO2 records')
% legend('\fontsize{12}La Jolla Station','\fontsize{12}Mauna Loa Station', '\fontsize{12}South Pole Observatory');
% 

%%% compute coherence

% this still seems to produce coherence of 1 everywhere
% cxy = 2*conj(MLOco2_ft(1:M+1,:)).*SPOco2_ft(1:M+1,:)/segment_length;
% cxy = cxy(2:length(cxy),:);
% C = abs(cxy)./sqrt(MLOco2_amp.*SPOco2_amp); 
% plot(C);
% 
% phase_C = atan2(-imag(cxy),real(cxy));

