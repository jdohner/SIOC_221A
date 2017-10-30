% file SIOC 221A HW 5
% 
% author Julia Dohner
%
% due date November 2, 2017

clear all; close all;

%% Compute two spectra

% generating 10,000 element dataset with Gaussian white noise
a=randn(10000,1);
% b(1)=b(1)=a(1); % this throws the error "The expression to the left of
% the equals sign is not a valid target for an assignment."
b(1) = a(1);
% generating a second dataset using autoregressive process
for i=2:length(a)
    b(i)=.5*b(i-1)+a(i);
end

% compute two spectra for white noise and autorecessive datasets by
% breaking the data up into segments

% reshaping data vector into 500x20 matrix
a_reshape = reshape(a,[500,20]);
b_reshape = reshape(b,[500,20]);
%a_test = a_reshape(:,1);

% % Question: when you reshape, does it calculate individual fft of each
% column? - Yes
f_a = fft(a_reshape);
f_b = fft(b_reshape);
% f_test = fft(a_test);
% f_firstcola = f_a(:,1);

% more rigorous way of ensuring calculating ffts of individual columns
% for i = 1:20
%     fft_a(:,i) = fft(a_reshape(:,i));
%     fft_b(:,i) = fft(b_reshape(:,i));
% end

% plotting spectra of a and b
N = 500; % length of each chunk of data
% Question: is this fine for the chunks in columns? probably no because
% amp_a becomes a 1x251 double
% amp_a = abs(f_a(1:N/2+1)).^2; % for even N
% amp_b = abs(f_b(1:N/2+1)).^2; % for even N

% pre-allocating for speed
amp_a = zeros(251,20);
amp_b = zeros(251,20);
for i = 1:20
    % Question: the +1 seems problematic because then I end up with 251x20
    % vectors
    amp_a(:,i) = abs(f_a(1:N/2+1)).^2; % for even N
    amp_b(:,i) = abs(f_b(1:N/2+1)).^2; % for even N
end


% summing over all spectra (taken from lecture 7 notes- seems weird)

% turning into row vectors to plot
amp_a = amp_a(:);
amp_b = amp_b(:);

% Question: how to create frequency vector?

% figure
% what is happening here is that each spectrum is individually calculated-
% they all look red. Then when they get plotted together like this they're
% just stuck next to each other, which makes it look like it's oscillating
semilogy(1:5020,amp_a, '-m', 1:5020, amp_b,'-b')
title('\fontsize{14}Spectra for white noise, autoregressive data');
xlabel('\fontsize{12}frequency');
ylabel('\fontsize{12}(units)^2 / frequency');
