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
valsSPO = textscan(dataSPO, '%f %f %f', ...
    'delimiter','\t \t');

fclose(dataMLO);
fclose(dataLJO);
fclose(dataSPO);

% format of .txt files is year, flag, co2 value
LJOyear = valsLJO{1};
LJOco2 = valsLJO{2};

MLOyear = valsMLO{1};
MLOco2 = valsMLO{2};

SPOyear = valsSPO{1};
SPOflag = valsSPO{2};
SPOco2 = valsSPO{3};

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
    if SPOflag(i) ~= 0
        SPOco2(i) = nan;
    end
end

% remove nan's
% TODO: if time, go back and change this to a linear interpolation
addpath('/Users/juliadohner/Documents/MATLAB/SIOC_221A/HW9/Inpaint_nans/Inpaint_nans');
MLOco2 = inpaint_nans(MLOco2);
LJOco2 = inpaint_nans(LJOco2);
SPOco2 = inpaint_nans(SPOco2);

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


% maybe time to segment!
% want to resolve an annual cycle

% start by taking segments of data of same length
startYear = LJOyear(1,1); %1.969041100000000e+03;
endYear =SPOyear(length(SPOyear),1);%1990.9562;
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

% % want to produce 3 segments with 50% overlap
% dataLength = length(LJOco2_2);
% M = 2; % number of segments splitting data into
% N = dataLength/M; % length of each chunk of data (aka segment length)
% 
% 
% LJOco2_3 = reshape(LJOco2_2,N,M);
% LJOco2_4 = LJOco2_2(N/2+1:length(LJOco2_2)-N/2);
% LJOco2_5 = [LJOco2_3 LJOco2_4];
% LJOco2_5 = detrend(LJOco2_5).*(hann(N));
% LJOco2_5 = fft(LJOco2_5);
% 
% MLOco2_3 = reshape(MLOco2_2,N,M);
% MLOco2_4 = MLOco2_2(N/2+1:length(MLOco2_2)-N/2);
% MLOco2_5 = [MLOco2_3 MLOco2_4];
% MLOco2_5 = detrend(MLOco2_5).*(hann(N));
% MLOco2_5 = fft(MLOco2_5);
% 
% % next step: calculate amplitude
% LJOco2_amp = (abs(LJOco2_5(1:N/2+1,:)).^2)/N;
% LJOco2_amp = LJOco2_amp(2:length(LJOco2_amp),:); % dump the mean
% LJOco2_mean = mean(LJOco2_amp,2); % take the mean across rows
% 
% MLOco2_amp = (abs(MLOco2_5(1:N/2+1,:)).^2)/N;
% MLOco2_amp = MLOco2_amp(2:length(MLOco2_amp),:); % dump the mean
% MLOco2_mean = mean(MLOco2_amp,2); % take the mean across rows
N = length(LJOco2_2);

LJOco2_5 = fft(detrend(LJOco2_2));
LJOco2_amp = (abs(LJOco2_5(1:N/2+1,:)).^2)/N;
LJOco2_amp = LJOco2_amp(2:length(LJOco2_amp),:);

MLOco2_5 = fft(detrend(MLOco2_2));
MLOco2_amp = (abs(MLOco2_5(1:N/2+1,:)).^2)/N;
MLOco2_amp = MLOco2_amp(2:length(MLOco2_amp),:); % dump the mean

SPOco2_5 = fft(detrend(SPOco2_2));
SPOco2_amp = (abs(SPOco2_5(1:N/2+1,:)).^2)/N;
SPOco2_amp = SPOco2_amp(2:length(SPOco2_amp),:); % dump the mean

frequency=(1:N/2)/(N/12);

% see harmonics on the annual cycle
figure
loglog(frequency,LJOco2_amp, '-r', frequency, MLOco2_amp, '-b', frequency, SPOco2_amp, '-g')
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}units')
title('\fontsize{16}Spectrum of Mauna Loa and La Jolla CO2 records')
legend('\fontsize{12}La Jolla Station','\fontsize{12}Mauna Loa Station', '\fontsize{12}South Pole Observatory');

%% compute coherence

