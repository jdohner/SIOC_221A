% file SIOC 221A HW 8
% 
% author Julia Dohner, with help from Annie Adelson
%
% due date November 30, 2017


%% aliasing

% fast sampling
f_sampling=1/(0.99349); % cycles per day
f_M2=24/12.42; % cycles per day
f_Nyquist=f_sampling/2;
M=floor(f_M2/f_Nyquist) % compute the integer ratio of the two frequencies.
alias = f_M2 - floor(f_M2/f_Nyquist)*f_Nyquist;
alias_period = 1/alias;
% Note: if M is odd then reset
if(rem(M,2)~=0) 
    alias=f_Nyquist-alias; 
end
%This calculation shows that M is 3 and the tidal energy aliases into a
% 2.3667 day period (alias_period)

% science sampling
f_sampling_sci=1/(20.86455); % cycles per day
f_Nyquist_sci=f_sampling_sci/2;
M_sci=floor(f_M2/f_Nyquist_sci) % compute the integer ratio of the two frequencies.
alias_sci = f_M2 - floor(f_M2/f_Nyquist_sci)*f_Nyquist_sci;
alias_period_sci = 1/alias_sci;
% Note: if M is odd then reset
if(rem(M,2)~=0) 
    alias_sci=f_Nyquist_sci-alias_sci; 
end
%This calculation shows that M is 80 and the tidal energy aliases into a
% 65.6178 day period (alias_period)