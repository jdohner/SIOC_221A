% file SIOC 221A HW 5
% 
% author Julia Dohner
%
% due date November 2, 2017

clear all; close all;

%% Compute two spectra - attempt via Lecture 7 Notes

% compute two spectra for white noise and autorecessive datasets by
% breaking the data up into segments

N = 10000; % length of each chunk of data
M = 20; % number of segments splitting data into
p = N/M; % datapoints per segment

% generating 10,000 element dataset with Gaussian white noise
a=randn(N,1);
b(1) = a(1);
% generating a second dataset using autoregressive process
for i=2:length(a)
    b(i)=.5*b(i-1)+a(i);
end

a_reshape = reshape(a,N/M,M);
b_reshape = reshape(b,N/M,M);
a_f = fft(a_reshape);
b_f = fft(b_reshape);

a_amp = abs(a_f(1:p/2+1,:)).^2; % aplitude of a = abs(f_t(from 1 to 500/2+1, all columns)^2
a_amp(2:p/2,:) = 2*a_amp(2:p/2,:)/p; % normalizing values in amp_a (except for 0)
% 
b_amp=abs(b_f(1:p/2+1,:)).^2; 
b_amp(2:p/2,:)=2*b_amp(2:p/2,:)/p;
a_amp = a_amp(2:251,:); % dumping the mean
a_amp = a_amp(1:249,:); % dumping last value
b_amp = b_amp(2:251,:); % dumping the mean
b_amp = b_amp(1:249,:); % dumping last value

frequency=(0:p/2)/p; % for N even
a_amp_mean = mean(a_amp,2);
b_amp_mean = mean(b_amp,2);

figure
semilogy(frequency,a_amp_mean, '-r', frequency, b_amp_mean, '-b')
xlabel('\fontsize{14}frequency')
ylabel('\fontsize{14}units')
title('\fontsize{16}Fourier Transform of Gaussian White Noise and Autoregressive Datasets')
legend('\fontsize{12}Gaussian white noise','Autoregressive data');


%% add error bars to your spectra

nu = 2*floor(N/p);
err_high = nu/chi2inv(0.05/2,nu);
err_low = nu/chi2inv(1-0.05/2,nu);
ratio_chi2 = err_high/err_low;

figure
% semilogy(frequency,mean_amp_a, '-r', frequency, mean_amp_b, '-b')
semilogy(0:p/2, a_amp_mean, '-b', [p/4 p/4],[err_low err_high]*a_amp_mean(p/4), '-r')
xlabel('\fontsize{12}2015')
ylabel('\fontsize{12}Pressure (dbar)')
title('\fontsize{14}2015 Scripps Pier Pressure')

% 0:M/2 is 10x1
% mean_amp_a, mean_amp_b are 251x1
% [M/4 M/4] is [5,5] (1x2)

% really confusing in end of lecture 7 notes because the definitions of M
% and p get switched
% last term is [1.4192, 1.1909] (1x2)

% seems to be ok, only plotting error bar halfway through our half
% spectrum?
figure
semilogy(0:p/2, b_amp_mean, '-b', [p/4 p/4],[err_low err_high]*b_amp_mean(p/4), '-r')
xlabel('\fontsize{12}2015')
ylabel('\fontsize{12}Pressure (dbar)')
title('\fontsize{14}2015 Scripps Pier Pressure')

%% monte carlo simulation

% bigMatrix is a matrix (500x4,000) containing 200 500x20 matrices
bigMatrix = randn(500,4000);
bigMatrix = fft(bigMatrix); % compute fft of each matrix

bigMatrix_amp=(abs(bigMatrix(1:p/2+1,:)).^2)/p; %needed to move the /p to this line from line below
bigMatrix_amp(2:p/2,:)=2*bigMatrix_amp(2:p/2,:);
bigMatrix_amp = bigMatrix_amp(2:251,:); % dumping the mean


for i=1:20:4000 %4000-19
    bigMatrix_mean_amp(:,floor(i/20)+1) = mean(bigMatrix_amp(:,i:i+19),2);
end

%turn mean amplitude result into a column vector
bigMatrix_pdf = bigMatrix_mean_amp(:);
bigMatrix_pdf = sort(bigMatrix_pdf,'descend'); 
err_low_data = prctile(bigMatrix_pdf,0.025);
err_high_data = prctile(bigMatrix_pdf,0.975);
ratio_monteCarlo = err_high_data/err_low_data;

bigMatrix_pdf_hist = histogram(bigMatrix_pdf); %histogram just a way of visualizing, just need info from raw values (have it once sorted)


%% repeat above steps, but using a Hanning window
% hanning: one dimensional, multiply by ones to make 2d
% should get roughly the same ratio when do hanning window

%so you multiply the hanning window before you fft the data. 
%so if you fft(X) and X has size NxM then you would fft(Y) where 
%Y = hanning(N)*(ones(1,M)).*X

% Hanning window treatment of bigMatrix in Question #4
bigMatrix_4 = randn(500,4000);
dataHan = fft(detrend(bigMatrix_4).*(hann(500)*ones(1,4000)));

dataHan_amp=(abs(dataHan(1:p/2+1,:)).^2)/p; %needed to move the /p to this line from line below
dataHan_amp(2:p/2,:)=2*dataHan_amp(2:p/2,:);
dataHan_amp = dataHan_amp(2:251,:); % dumping the mean

for i=1:20:4000 %4000-19
    dataHan_mean_amp(:,floor(i/20)+1) = mean(dataHan_amp(:,i:i+19),2);
end

%turn mean amplitude result into a column vector
dataHan_pdf = dataHan_mean_amp(:);
dataHan_pdf = sort(dataHan_pdf,'descend'); 
err_low_data = prctile(dataHan_pdf,0.025);
err_high_data = prctile(dataHan_pdf,0.975);
ratio_Hann = err_high_data/err_low_data;

dataHan_pdf_hist = histogram(dataHan_pdf); 