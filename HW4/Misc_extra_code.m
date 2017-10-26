%% THE CODE THAT WORKS

% time differences
time = double(time); 
t_diff = diff(time);
%t_diff2 = t_diff/60;
big_diff = max(t_diff);


%kAugustStart = floor((7/12)*82237); % =47971
cutoff = find(t_diff(1:82236)>t_diff(1),1);
pressureAugust = pressure(1:cutoff-1);
timeAugust = time(1:cutoff-1);
date0=datenum(1970,1,1); % give reference date (first date)
%timeAugust = timeAugust/86400+date0;

%% CODE that I want to make work:

kAugustStart = floor((7/12)*82237); % =47971

% time differences
time = double(time); 
t_diff = diff(time);
%t_diff2 = t_diff/60;
%big_diff = max(t_diff);

t_diff_pastAugust = t_diff(kAugustStart:length(t_diff)); % new vectors of just t_diff beyond august
big_diff_pastAugust = max(t_diff_pastAugust);


% bump the cutoff to be in the summer ~august

cutoff = find(t_diff_pastAugust>t_diff_pastAugust(1),1);
cutoff_full = cutoff+kAugustStart;
pressureAugust = pressure(kAugustStart:cutoff_full-1);
timeAugust = time(kAugustStart:cutoff_full-1);
date0=datenum(1970,1,1); % give reference date (first date)
timeAugust = timeAugust/86400+date0;

%% old old code that wasn't working


% %% Pier pressure record for a summer month of 2015
% 
% % 30 days roughly near August 2015 
% 
% % Starting 7/12 of the way through the time record (82237 measurements).
% % If a measurement is taken every 361 seconds, then 30 days into the record
% % should be roughly the next 7180 measurements (2592000/361). My record is 
% % 30 days long.
% date0=datenum(1970,1,1); % give reference date (first date)
% kAugustStart = floor((7/12)*82237); 
% kAugustEnd = kAugustStart + 7180; 
% timeAugust = time(kAugustStart:kAugustEnd, 1);
% timeAugust = double(timeAugust)/86400+date0; % in units of days (conversion: seconds/day)
% %timeAugustHours = double(timeAugust)/3600 + date0;
% pressureAugust = pressure(kAugustStart:kAugustEnd, 1);
% % 7181 points in record at 361 second intervals 
% % total duration is 7181*361/(24*3600) days
% 
% %  %% time differences
% % %timeCheck = double(time);
% % t_diff = diff(timeAugust);
% % t_diff2 = t_diff/60;
% % big_diff = max(t_diff2);