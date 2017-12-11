% file SIOC 221A HW 9
% 
% author Julia Dohner
%
% due date December 13, 2017

% analysis of daily flask-collected co2 data at Mauna Loa and Point Barrow
% stations
% from 1961-1963.74, 1965.45-1967.7 (avoiding gaps in PTB data)
%

clear all; close all;

%% load CO2 data - daily

dataMLO = fopen('daily_data/daily_in_situ_co2_mlo_JLD.txt');
dataPTB = fopen('daily_data/daily_in_situ_co2_ptb_JLD.txt');

valsMLO = textscan(dataMLO, '%f %f', ...
    'delimiter','\t');
valsPTB = textscan(dataPTB, '%f %f %f', ...
    'delimiter','\t');

fclose(dataMLO);
fclose(dataPTB);


% format of .txt files is year, co2 value
MLOyear = valsMLO{1};
MLOco2 = valsMLO{2};

PTByear = valsPTB{1};
PTBflag = valsPTB{2};
PTBco2 = valsPTB{3};

% remove flagged data
for i = 1:length(PTBco2)
    if PTBflag(i) ~= 0
        PTBco2(i) = nan;
    end
end


% remove nan's
addpath('/Users/juliadohner/Documents/MATLAB/SIOC_221A/HW9/Inpaint_nans/Inpaint_nans');
MLOyear = inpaint_nans(MLOyear);
MLOco2 = inpaint_nans(MLOco2);
PTByear = inpaint_nans(PTByear);
PTBco2 = inpaint_nans(PTBco2);



%% inspect spacing of data

% inspect time spacing between measurements
MLO_t_diff = diff(MLOyear);
PTB_t_diff = diff(PTByear);

% figure
% plot(1:length(MLOyear)-1,MLO_t_diff);
% figure
% plot(1:length(LJOyear)-1,LJO_t_diff);
% figure
% plot(1:length(SPOyear)-1,SPO_t_diff);

minDiffPTB = min(PTB_t_diff);
maxDiffPTB = max(PTB_t_diff);
minDiffMLO = min(MLO_t_diff);
maxDiffMLO = max(MLO_t_diff);

% too big differences in time increments of daily LJO flask sampling, will try 
% monthly data instead

% update: time differences range between 0.0767 and 0.0850 years for both 
% MLO and LJO records (difference of 0.0083) or ~3 days, deemed small 
% enough difference in time increments
% read: counts in my book as even spacing

%% shorten data to same lengths

% start records at earliest date in 1968 (LJO's latest start year)
% startYear = LJOyear(1,1); % 1.969041100000000e+03 (latest start)
% endYear = SPOyear(length(SPOyear),1); % earliest end
% startIndex_LJO = find(LJOyear == startYear);
% startIndex_MLO = find(MLOyear == startYear);
% startIndex_SPO = find(SPOyear == startYear);
% endIndex_LJO = find(LJOyear == endYear);
% endIndex_MLO = find(MLOyear == endYear);
% endIndex_SPO = find(SPOyear == endYear);

% from 1961-1963.74, 1965.45-1967.7 (avoiding gaps in PTB data)
startYear1 = 1962.12; % looks like resolution improves in LJO after 1985
endYear1 = 1963.74;
startIndex1_PTB = find(floor(PTByear*100) == startYear1*100, 1);
startIndex1_MLO = find(floor(MLOyear*100) == startYear1*100,1);
% endIndex1_PTB = find(floor(PTByear*100) == endYear1*100,1);
% endIndex1_MLO = find(floor(MLOyear*100) == endYear1*100,1);

startYear2 = 1966.49;
endYear2 = 1967.70;
startIndex2_PTB = find(floor(PTByear*100) == startYear2*100, 1);
startIndex2_MLO = find(floor(MLOyear*100) == startYear2*100,1);
endIndex2_PTB = find(floor(PTByear*100) == endYear2*100,1);
endIndex2_MLO = find(floor(MLOyear*100) == endYear2*100,1);

% some data missing, so interpolate to get vectors of same length

dt = 1/365; % daily increment in unit of years

year1 = startYear1:dt:endYear1;
year2 = startYear2:dt:endYear2;
PTBco2_1 = interp1(PTByear,PTBco2,year1,'spline');
MLOco2_1 = interp1(MLOyear,MLOco2,year1,'spline');

PTBco2_2 = interp1(PTByear,PTBco2,year2,'spline');
MLOco2_2 = interp1(MLOyear,MLOco2,year2,'spline');




% PTBco2_firstHalf = PTBco2(startIndex1_PTB:endIndex1_PTB);
% PTByear_firstHalf = PTByear(startIndex1_PTB:endIndex1_PTB);
% MLOco2_firstHalf = MLOco2(startIndex1_MLO:endIndex1_MLO);
% MLOyear_firstHalf = MLOyear(startIndex1_MLO:endIndex1_MLO);
% 
% PTBco2_secondHalf = PTBco2(startIndex2_PTB:endIndex2_PTB);
% PTByear_secondHalf = PTByear(startIndex2_PTB:endIndex2_PTB);
% MLOco2_secondHalf = MLOco2(startIndex2_MLO:endIndex2_MLO);
% MLOyear_secondHalf = MLOyear(startIndex2_MLO:endIndex2_MLO);
% 
% PTBco2_2 = vertcat(PTBco2_firstHalf,PTBco2_secondHalf);
% PTByear_2 = vertcat(PTByear_firstHalf,PTByear_secondHalf);
% MLOco2_2 = vertcat(MLOco2_firstHalf, MLOco2_secondHalf);
% MLOyear_2 = vertcat(MLOyear_firstHalf, MLOyear_secondHalf);

% plot timeseries
figure('name','Atmospheric CO2 Timeseries');
plot(year1,PTBco2_1, '-r', year1,MLOco2_1,'-b', year2,PTBco2_2,'-r',year2,MLOco2_2,'-b');
xlabel('\fontsize{14}year')
ylabel('\fontsize{14}ppm')
title('\fontsize{16}Atmospheric CO2 Records')
legend('\fontsize{12}Point Barrow','\fontsize{12}Mauna Loa', 'Location','northwest');


%% compute spectra (from in-class coherence example)

segment_length = 44; % length of each chunk of data (aka segment length)

% cut records down by two to make more easily divisible into segments
length_1 = floor(length(MLOco2_1)/segment_length)*segment_length;
length_2 = floor(length(MLOco2_2)/segment_length)*segment_length;

year1 = year1(1:length_1);
MLOco2_1 = MLOco2_1(1:length_1);
PTBco2_1 = PTBco2_1(1:length_1);

year2 = year2(1:length_2);
MLOco2_2 = MLOco2_2(1:length_2);
PTBco2_2 = PTBco2_2(1:length_2);

N1 = length_1;
Nseg1 = N1/segment_length; % number of segments splitting data into
N2 = length_2;
Nseg2 = N2/segment_length;
M = segment_length/2;

Nseg = Nseg1 + Nseg2;

% using three segments with 50% overlap (concatenating in lines below)
% LJO_use=[reshape(LJOco2_2,segment_length,Nseg) reshape(LJOco2_2(M+1:end-M),segment_length,Nseg-1)];
% MLO_use=[reshape(MLOco2_2,segment_length,Nseg) reshape(MLOco2_2(M+1:end-M),segment_length,Nseg-1)];
% SPO_use=[reshape(SPOco2_2,segment_length,Nseg) reshape(SPOco2_2(M+1:end-M),segment_length,Nseg-1)];

% changing to non-overlapping so can fiddle with number of segments:
PTB_use1=reshape(PTBco2_1,segment_length,Nseg1);
PTB_use2=reshape(PTBco2_2,segment_length,Nseg2);
MLO_use1=reshape(MLOco2_1,segment_length,Nseg1);
MLO_use2=reshape(MLOco2_2,segment_length,Nseg2);

% take fft of first and second bits of data
PTB_ft1=fft(detrend(PTB_use1).*(hann(segment_length)*ones(1,Nseg1)));
PTB_ft2=fft(detrend(PTB_use2).*(hann(segment_length)*ones(1,Nseg2)));

MLO_ft1=fft(detrend(MLO_use1).*(hann(segment_length)*ones(1,Nseg1)));
MLO_ft2=fft(detrend(MLO_use2).*(hann(segment_length)*ones(1,Nseg2)));

PTB_ft = [PTB_ft1 PTB_ft2];
MLO_ft = [MLO_ft1 MLO_ft2];

% question: mean is still in here. Do I want it?
PTB_spec=sum(abs(PTB_ft(1:M+1,:)).^2,2)/N1; % sum over all spectra (2nd dim)
PTB_spec(2:end)=PTB_spec(2:end)*2; % multiply by 2 to make up for lost energy

MLO_spec=sum(abs(MLO_ft(1:M+1,:)).^2,2)/N1;
MLO_spec(2:end)=MLO_spec(2:end)*2;


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
loglog(frequency,PTB_spec, '-r', [.2 .2],[err_low err_high]*PTB_spec(10), ...
    frequency, MLO_spec, '-b',[.1 .1],[err_low err_high]*MLO_spec(10))
xlabel('\fontsize{14}cycles per year')
ylabel('\fontsize{14}ppm^2/cpy')
title('\fontsize{16}Power Spectra of CO2 Records')
legend('\fontsize{12}La Jolla Station','\fontsize{12}Mauna Loa Station', 'Location','northeast');

%% compute coherence

% compute cross covariance of PTB/MLO
ccLM=sum(PTB_ft(1:M+1,:).*conj(MLO_ft(1:M+1,:)),2)/N1;
ccLM(2:end)=ccLM(2:end)*2;
% compute coherence
C_LM=abs(ccLM)./sqrt(PTB_spec.*MLO_spec);
phase_LM = atan2(-imag(ccLM),real(ccLM));
deltaPhase_LM = sqrt((1-C_LM.^2)./(abs(C_LM).^2*2*Nseg1)); 

% % compute cross covariance of MLO/SPO
% ccMS=sum(MLO_ft(1:M+1,:).*conj(SPO_ft(1:M+1,:)),2)/N1;
% ccMS(2:end)=ccMS(2:end)*2;
% % compute coherence
% C_MS=abs(ccMS)./sqrt(MLO_spec.*SPO_spec);
% phase_MS = atan2(-imag(ccMS),real(ccMS));
% deltaPhase_MS = sqrt((1-C_MS.^2)./(abs(C_MS).^2*2*Nseg1));

%% uncertainty estimate 

alpha = 0.05; % 95% confidence level
gamma_threshold= sqrt(1-alpha^(1/(Nseg1-1)));


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