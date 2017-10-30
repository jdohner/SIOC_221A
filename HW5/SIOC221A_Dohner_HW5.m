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
f_a = fft(a_reshape);
f_b = fft(b_reshape);
% for i = 1:floor(N/M);
%     amp_a(:,i) = abs(f_a(1:(N/2+1))).^2; % for even N
%     amp_b(:,i) = abs(f_b(1:(N/2+1))).^2; % for even N
% end

amp_a = abs(f_a(1:p/2+1,:)).^2; % aplitude of a = abs(f_t(from 1 to 500/2+1, all columns)^2
amp_a(2:p/2,:) = 2*amp_a(2:p/2,:)/p; % normalizing values in amp_a (except for 0)
% 
amp_b=abs(f_b(1:p/2+1,:)).^2; 
amp_b(2:p/2,:)=2*amp_b(2:p/2,:)/p;

frequency=(0:p/2)/p; % for N even
mean_amp_a = mean(amp_a,2);
mean_amp_b = mean(amp_b,2);
figure
semilogy(frequency,mean_amp_a, '-r', frequency, mean_amp_b, '-b')
xlabel('\fontsize{12}2015')
ylabel('\fontsize{12}Pressure (dbar)')
title('\fontsize{14}2015 Scripps Pier Pressure')

% %% Compute two spectra - my initial attempt
% 
% % generating 10,000 element dataset with Gaussian white noise
% a=randn(10000,1);
% b(1) = a(1);
% % generating a second dataset using autoregressive process
% for i=2:length(a)
%     b(i)=.5*b(i-1)+a(i);
% end
% 
% % compute two spectra for white noise and autorecessive datasets by
% % breaking the data up into segments
% 
% % reshaping data vector into 500x20 matrix
% a_reshape = reshape(a,[500,20]); % suggestion for getting red, white spectra
% b_reshape = reshape(b,[500,20]);
% 
% f_a = fft(a_reshape);
% f_b = fft(b_reshape);
% 
% % plotting spectra of a and b
% N = 500; % length of each chunk of data
% 
% % pre-allocating for speed
% amp_a = zeros(251,20);
% amp_b = zeros(251,20);
% for i = 1:20
%     amp_a(:,i) = abs(f_a(1:N/2+1)).^2; % for even N
%     amp_b(:,i) = abs(f_b(1:N/2+1)).^2; % for even N
% end
% 
% % take mean across rows
% mean_amp_a = mean(amp_a,2)/(N/2); 
% mean_amp_b = mean(amp_b,2)/(N/2); 
% % multiplying by 2 bc taking half frequencies, divide by N to normalize
% % technically don't do to freq 0, but won't be looking at freq 0
% 
% figure
% semilogy(0:250,mean_amp_a, 0:250, mean_amp_b, '-b');
% title('\fontsize{14}Spectra for white noise, autoregressive data');
% xlabel('\fontsize{12}frequency'); % freq is arbitrary, technically cycles per length of data
% ylabel('\fontsize{12}(units)^2 / frequency'); % y is arbitrary energy, raw data squared, divide by freq

% Helpful meeting notes:
%
% 2 in mean(x,2) is taking mean across second dimension (each row)
% default is dim 1, which is each column
% 250 freqs, 20 datasets
% can think of 20 thermistors on pier, ft each time series, average
% together to get typical spectrum for whole record

% seems like there's a lot of noise, go on to next section, look at error
% bars, see if wiggles are consistent with error bars, if not, look at code
%could also make dataset larger (200 instead of 20)

% summing: want energy in each segment (lect 7 notes)
% mean: want typical energy
% need to divide by length of segment to get the energy to work out
% (normalize appropriately)

%% add error bars to your spectra

nu = 2*floor(N/p);
err_low = nu/chi2inv(0.05/2,nu);
err_high = nu/chi2inv(1-0.05/2,nu);

figure
% semilogy(frequency,mean_amp_a, '-r', frequency, mean_amp_b, '-b')
semilogy(0:p/2, mean_amp_a, '-b', [p/4 p/4],[err_low err_high]*mean_amp_a(p/4), '-r')
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
semilogy(0:p/2, mean_amp_b, '-b', [p/4 p/4],[err_low err_high]*mean_amp_b(p/4), '-r')
xlabel('\fontsize{12}2015')
ylabel('\fontsize{12}Pressure (dbar)')
title('\fontsize{14}2015 Scripps Pier Pressure')

