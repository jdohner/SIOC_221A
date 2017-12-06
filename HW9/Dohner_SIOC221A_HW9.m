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


% maybe time to segment!
% want to resolve an annual cycle

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

%% create 3 segments with 50% overlap


N = length(LJOco2_2);
Nseg = 2; % number of segments splitting data into
segmentLength = N/Nseg; % length of each chunk of data (aka segment length)
M = segmentLength/2;

LJOco2_3 = reshape(LJOco2_2,segmentLength,Nseg);
LJOco2_4 = LJOco2_2(M+1:length(LJOco2_2)-M);
LJOco2_5 = [LJOco2_3 LJOco2_4];
LJOco2_5 = detrend(LJOco2_5).*(hann(segmentLength)*ones(1,3)); % 3 = number of columns
LJOco2_ft = fft(LJOco2_5);

MLOco2_3 = reshape(MLOco2_2,segmentLength,Nseg);
MLOco2_4 = MLOco2_2(M+1:length(MLOco2_2)-M);
MLOco2_5 = [MLOco2_3 MLOco2_4];
MLOco2_5 = detrend(MLOco2_5).*(hann(segmentLength)*ones(1,3));
MLOco2_ft = fft(MLOco2_5);

SPOco2_3 = reshape(SPOco2_2,segmentLength,Nseg);
SPOco2_4 = SPOco2_2(M+1:length(SPOco2_2)-M);
SPOco2_5 = [SPOco2_3 SPOco2_4];
SPOco2_5 = detrend(SPOco2_5).*(hann(segmentLength)*ones(1,3));
SPOco2_ft = fft(SPOco2_5);

% next step: calculate amplitude
LJOco2_amp = 2*(abs(LJOco2_ft(1:M+1,:)).^2)/segmentLength;
LJOco2_amp = LJOco2_amp(2:length(LJOco2_amp),:); % dump the mean
LJOco2_mean = mean(LJOco2_amp,2); % take the mean across rows

MLOco2_amp = 2*(abs(MLOco2_ft(1:M+1,:)).^2)/segmentLength;
MLOco2_amp = MLOco2_amp(2:length(MLOco2_amp),:); % dump the mean
MLOco2_mean = mean(MLOco2_amp,2); % take the mean across rows

SPOco2_amp = 2*(abs(SPOco2_ft(1:M+1,:)).^2)/segmentLength;
SPOco2_amp = SPOco2_amp(2:length(SPOco2_amp),:); % dump the mean
SPOco2_mean = mean(SPOco2_amp,2); % take the mean across rows

frequency=(1:M)/(segmentLength/12);
frequency = frequency';

% plot spectra
% uncertainty estimate
nu = 2*1; % DOF = 2*number of segments
% do I need to segment?
err_high_MLO = nu/chi2inv(0.05/2,nu);
err_low_MLO = nu/chi2inv(1-0.05/2,nu);
ratio_chi2_MLO = err_high_MLO/err_low_MLO;

% see harmonics on the annual cycle

figure('name','spectra with 3 segments, 50% overlap');
loglog(frequency,LJOco2_mean, '-r', frequency, MLOco2_mean, '-b',[.1 .1],[err_low_MLO err_high_MLO]*MLOco2_amp(100),frequency, SPOco2_mean, '-g')
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}ppm^2/cpy')
title('\fontsize{16}Spectrum of Mauna Loa and La Jolla CO2 records')
legend('\fontsize{12}La Jolla Station','\fontsize{12}Mauna Loa Station', '\fontsize{12}South Pole Observatory');


%% compute coherence

cxy = 2*conj(MLOco2_ft(1:M+1,:)).*SPOco2_ft(1:M+1,:)/segmentLength;
cxy = cxy(2:length(cxy),:);
C = abs(cxy)./sqrt(MLOco2_amp.*SPOco2_amp); % comparing new method to old
plot(C);
% cxy = cxy(2:length(cxy),:);
% C=abs(cxy)./sqrt(MLOco2_amp.*SPOco2_amp);
% plot(C)

% compute cross-spectrum
% cxy=sum(MLOco2_ft(2:M+1,:).*conj(SPOco2_ft(2:M+1,:)),2)/segmentLength;
% cxy(2:end)=cxy(2:end)*2; % multiply cospectrum by 2 (bc multiplied spec *2)
% C=abs(cxy)./sqrt(MLOco2_mean.*SPOco2_mean);
% phase_C = atan2(-imag(cxy),real(cxy));
% figure
% plot(C);

% segment_length=500;
% N=length(x);
% M=segment_length/2; % define this value - 250
% Nseg=N/segment_length; % number of segments

% %% compute spectra
% 
% N = length(LJOco2_2);
% 
% LJOco2_ft = fft(detrend(LJOco2_2));
% LJOco2_amp = (abs(LJOco2_ft(1:N/2+1,:)).^2)/N;
% LJOco2_amp = LJOco2_amp(2:length(LJOco2_amp),:);
% 
% MLOco2_ft = fft(detrend(MLOco2_2));
% MLOco2_amp = 2*(abs(MLOco2_ft(1:N/2+1,:)).^2)/N;
% MLOco2_amp = MLOco2_amp(2:length(MLOco2_amp),:); % dump the mean
% 
% SPOco2_ft = fft(detrend(SPOco2_2));
% SPOco2_amp = 2*(abs(SPOco2_ft(1:N/2+1,:)).^2)/N;
% SPOco2_amp = SPOco2_amp(2:length(SPOco2_amp),:); % dump the mean

% frequency = frequency';

%cxy=conj(MLOco2_ft(2:end/2)).*SPOco2_ft(2:end/2)/N;


%% compute coherence

% cxy=conj(MLOco2_amp(1:end/2)).*SPOco2_amp(1:end/2); 
% cxy(2:end)=2*cxy(2:end);
% C=abs(cxy)./sqrt(MLOco2_amp.*SPOco2_amp); % matrix dims must agree
% plot(C)


% 
% % cab is covariance between a and b
% alpha = .05;
% nd=10; % # of segments
% p=1-alpha;
% delta_phase2 = sqrt((1-cab.?2)./(abs(cab).?2*2*nd));

% should I segment?
% making sense of axes in coherence
% annie's ADCP: depth integrate or at one depth? - can do either
% should I expect coherence to be 1 everywhere because I didn't segment?

% TODO: segment my monthly data
% try this out for daily measurements to see if can get time lag