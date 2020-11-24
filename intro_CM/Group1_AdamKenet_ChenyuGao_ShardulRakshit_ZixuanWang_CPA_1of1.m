% Coding Project Assignment
% 1 of 1

% Group 1: Adam Kenet, Chenyu Gao, Shardul Rakshit, Zixuan Wang
% Introduction to Computational Medicine: The Physiome
% November 18, 2020

% NOTE: Make sure all files are in same directory or paths are correct.
%       '*_ABP.txt'
%       '*n.txt' NOTE: manually remove row of units (second row) in file before importing
%       wabp.m
%       abpfeature.m
%       jSQI.m
%       estimateCO_v3.m
%       'est#_*.m'
%       estimate_parlikar.m

%% A: (Q1) Patient 20: First 20 ABP beats after 10 hrs
clear variables; clc;

[ABP,onsets,EoS_m1,EoS_m2] = abpANDfeatures('s00020-2567-03-30-17-47_ABP.txt',10,20);

% plot
figure(1)
plot(ABP) % ABP waveform
hold on;
p1=plot(onsets(:,1), onsets(:,2), '*') % onset
p2=plot(EoS_m1(:,1), EoS_m1(:,2), 'X') % end of systole (method 1)
p3=plot(EoS_m2(:,1), EoS_m2(:,2), 'O') % end of systole (method 2)
xlabel('Sample Number (125 = 1 second)', 'FontSize',12);
ylabel('ABP (mmHg)', 'FontSize',12);
title('Patient 20: First 20 ABP Pulses after 10 hours', 'FontSize',12);
legend([p1,p2,p3],'onset point','end of systole (method 1)','end of systole (method 2)');

%% B: (Q1) Patient 20: First 20 ABP beats after 11 hrs
clear variables; clc;

[ABP,onsets,EoS_m1,EoS_m2] = abpANDfeatures('s00020-2567-03-30-17-47_ABP.txt',11,20);

% plot
figure(2)
plot(ABP) % ABP waveform
hold on;
p1=plot(onsets(:,1), onsets(:,2), '*') % onset
p2=plot(EoS_m1(:,1), EoS_m1(:,2), 'X') % end of systole (method 1)
p3=plot(EoS_m2(:,1), EoS_m2(:,2), 'O') % end of systole (method 2)
xlabel('Sample Number (125 = 1 second)', 'FontSize',12);
ylabel('ABP (mmHg)', 'FontSize',12);
title('Patient 20: First 20 ABP Pulses after 11 hours', 'FontSize',12);
legend([p1,p2,p3],'onset point','end of systole (method 1)','end of systole (method 2)');

%% C: (Q2) Patient 214: First 20 ABP beats after 18 hrs
clear variables; clc;

[ABP,onsets,EoS_m1,EoS_m2] = abpANDfeatures('s00214-3084-11-28-16-23_ABP.txt',18,20);

% plot
figure(3)
plot(ABP) % ABP waveform
hold on;
p1=plot(onsets(:,1), onsets(:,2), '*') % onset
p2=plot(EoS_m1(:,1), EoS_m1(:,2), 'X') % end of systole (method 1)
p3=plot(EoS_m2(:,1), EoS_m2(:,2), 'O') % end of systole (method 2)
xlabel('Sample Number (125 = 1 second)', 'FontSize',12);
ylabel('ABP (mmHg)', 'FontSize',12);
title('Patient 214: First 20 ABP Pulses after 18 hours', 'FontSize',12);
legend([p1,p2,p3],'onset point','end of systole (method 1)','end of systole (method 2)');

%% D: (Q2) Patient 5237: First 20 ABP beats after 15 hrs
clear variables; clc;

[ABP,onsets,EoS_m1,EoS_m2] = abpANDfeatures('s05237-2926-10-24-11-50_ABP.txt',15,20);

% plot
figure(4)
plot(ABP) % ABP waveform
hold on;
p1=plot(onsets(:,1), onsets(:,2), '*') % onset
p2=plot(EoS_m1(:,1), EoS_m1(:,2), 'X') % end of systole (method 1)
p3=plot(EoS_m2(:,1), EoS_m2(:,2), 'O') % end of systole (method 2)
xlabel('Sample Number (125 = 1 second)', 'FontSize',12);
ylabel('ABP (mmHg)', 'FontSize',12);
title('Patient 5237: First 20 ABP Pulses after 15 hours', 'FontSize',12);
legend([p1,p2,p3],'onset point','end of systole (method 1)','end of systole (method 2)');

%% E: (Q3) Patient 20: Liljestrand Estimation and C2 Calibration of CO, PP, MAP, HR
clear variables; clc;

num_filename = 's00020-2567-03-30-17-47n.txt';
ABP_filename = 's00020-2567-03-30-17-47_ABP.txt';
estID = 5; % Liljestrand estimator
[feature_times,features,valsATtco,~] = estimateANDcalibrate(num_filename,ABP_filename,estID);

% get features
CO = features(:,1); % cardiac output [Lpm]
PP = features(:,2); % pulse pressure [mmHg]
MAP = features(:,3); % mean arterial pressure [mmHg]
HR = features(:,4); % heart rate [beats/min]

% get values at experimental TCO
TCO_times = valsATtco(:,1); % [hr]
TCO_values = valsATtco(:,2); % [Lpm]
PP_values = valsATtco(:,3); % [mmHg]
MAP_values = valsATtco(:,4); % [mmHg]
HR_values = valsATtco(:,5); % [beats/min]

% plot CO
figure(5)
subplot(4,1,1)
plot(feature_times,CO) % C2 calibrated CO [Lpm]
hold on;
stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
ylabel('CO', 'FontSize',12);

% plot PP
subplot(4,1,2)
plot(feature_times,PP) % pulse pressure [mmHg]
hold on;
stem(TCO_times,PP_values,'filled') % PP values [mmHg]
ylabel('PP', 'FontSize',12);

% plot MAP
subplot(4,1,3)
plot(feature_times,MAP) % mean arterial pressure [mmHg]
hold on;
stem(TCO_times,MAP_values,'filled') % MAP values [mmHg]
ylabel('MAP', 'FontSize',12);

% plot HR
subplot(4,1,4)
plot(feature_times,HR) % heart rate [beats/min]
hold on;
stem(TCO_times,HR_values,'filled') % HR values [beats/min]
ylabel('HR', 'FontSize',12);
xlabel('time (hours)', 'FontSize',12);
sgtitle('Patient 20: Liljestrand Estimation and C2 Calibration')

%% F: (Q4/5) Patient 20: Comparison of Different CO Estimation Algorithms
clear variables; clc;

num_filename = 's00020-2567-03-30-17-47n.txt';
ABP_filename = 's00020-2567-03-30-17-47_ABP.txt';

estIDs = [1,2,3,4,5,6,7,10,2007];
estNames = ["MAP","Windkessel","SA","SA Warner","Liljestrand","Herd","SA Wesseling","RC Decay","Parlikar"];

figure(6)
set(gcf, 'Position', [100, 100, 1300, 700])
for i = 1 : length(estIDs)
    estID = estIDs(i);
    [feature_times,features_estimated,valsATtco,~] = estimateANDcalibrate(num_filename,ABP_filename,estID);
    
    % C2 calibrate
    CO_estimated = features_estimated(:,1); % C2 calibrated cardiac output [Lpm]
    
    % RMSNE
    TCO_times = valsATtco(:,1); % times of experimental TCO [hr]
    TCO_values = valsATtco(:,2); % values of experimental TCO [Lpm]
    rmsne = RMSNE(TCO_times, TCO_values, feature_times, CO_estimated); % calculate RMSNE
    all_rmsne(i) = rmsne; % each column (=estID) is RMSNE of different estimator
    
    % plot
    subplot(5,2,i)
    plot(feature_times,CO_estimated) % CO from MAP estimator [Lpm]
    hold on;
    stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
    ylim([0 15]);
    title(estNames(i), 'FontSize',12);
end

% Robust Parlikar
[feature_times,features,valsATtco] = RobustParlikar(num_filename,ABP_filename);
% RMSNE
CO_estimated = features(:,1); % cardiac output [Lpm]
TCO_times = valsATtco(:,1); % times of experimental TCO [hr]
TCO_values = valsATtco(:,2); % values of experimental TCO [Lpm]
rmsne_RP = RMSNE(TCO_times, TCO_values, feature_times, CO_estimated); % calculate RMSNE

subplot(5,2,10)
plot(feature_times,CO_estimated) % CO from MAP estimator [Lpm]
hold on;
stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
ylim([0 15]);
xlabel('Time (hr)', 'FontSize',12)
title('Robust Parlikar', 'FontSize',12);

% fix axes labels
subplot(5,2,1); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,3); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,5); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,7); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,9); ylabel('CO (Lpm)', 'FontSize',12); xlabel('Time (hr)', 'FontSize',12);
sgtitle('Patient 20: Comparison of Different CO Estimation Algorithms')

figure(7)
bar([all_rmsne(1), all_rmsne(2), all_rmsne(3), all_rmsne(4), all_rmsne(5), all_rmsne(6), all_rmsne(7), all_rmsne(8), all_rmsne(9), rmsne_RP])
ylabel('RMSNE', 'FontSize',12);
xticklabels({'MAP','Windkessel','SA','SA Warner','Liljestrand','Herd','SA Wesseling', 'RC Decay', 'Parlikar', 'Robust Parlikar'})
xtickangle(45)
title('Patient 20: RMSNE Different CO Estimation Algorithms', 'FontSize',12);

%% G: (Q4/5) Patient 214: Comparison of Different CO Estimation Algorithms
clear variables; clc;

num_filename = 's00214-3084-11-28-16-23n.txt';
ABP_filename = 's00214-3084-11-28-16-23_ABP.txt';

estIDs = [1,2,3,4,5,6,7,10,2007];
estNames = ["MAP","Windkessel","SA","SA Warner","Liljestrand","Herd","SA Wesseling","RC Decay","Parlikar"];

figure(8)
set(gcf, 'Position', [100, 100, 1300, 700])
for i = 1 : length(estIDs)
    estID = estIDs(i);
    [feature_times,features_estimated,valsATtco,~] = estimateANDcalibrate(num_filename,ABP_filename,estID);
    
    % C2 calibrate
    CO_estimated = features_estimated(:,1); % C2 calibrated cardiac output [Lpm]
    
    % RMSNE
    TCO_times = valsATtco(:,1); % times of experimental TCO [hr]
    TCO_values = valsATtco(:,2); % values of experimental TCO [Lpm]
    rmsne = RMSNE(TCO_times, TCO_values, feature_times, CO_estimated); % calculate RMSNE
    all_rmsne(i) = rmsne; % each column (=estID) is RMSNE of different estimator
    
    % plot
    subplot(5,2,i)
    plot(feature_times,CO_estimated) % CO from MAP estimator [Lpm]
    hold on;
    stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
    ylim([0 15]);
    title(estNames(i), 'FontSize',12);
end

% Robust Parlikar
[feature_times,features,valsATtco] = RobustParlikar(num_filename,ABP_filename);
% RMSNE
CO_estimated = features(:,1); % cardiac output [Lpm]
TCO_times = valsATtco(:,1); % times of experimental TCO [hr]
TCO_values = valsATtco(:,2); % values of experimental TCO [Lpm]
rmsne_RP = RMSNE(TCO_times, TCO_values, feature_times, CO_estimated); % calculate RMSNE

subplot(5,2,10)
plot(feature_times,CO_estimated) % CO from MAP estimator [Lpm]
hold on;
stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
ylim([0 15]);
xlabel('Time (hr)', 'FontSize',12)
title('Robust Parlikar', 'FontSize',12);

% fix axes labels
subplot(5,2,1); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,3); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,5); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,7); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,9); ylabel('CO (Lpm)', 'FontSize',12); xlabel('Time (hr)', 'FontSize',12);
sgtitle('Patient 214: Comparison of Different CO Estimation Algorithms')

figure(9)
bar([all_rmsne(1), all_rmsne(2), all_rmsne(3), all_rmsne(4), all_rmsne(5), all_rmsne(6), all_rmsne(7), all_rmsne(8), all_rmsne(9), rmsne_RP])
ylabel('RMSNE', 'FontSize',12);
xticklabels({'MAP','Windkessel','SA','SA Warner','Liljestrand','Herd','SA Wesseling', 'RC Decay', 'Parlikar', 'Robust Parlikar'})
xtickangle(45)
title('Patient 214: RMSNE Different CO Estimation Algorithms', 'FontSize',12);

%% H: (Q4/5) Patient 5237: Comparison of Different CO Estimation Algorithms
clear variables; clc;

num_filename = 's05237-2926-10-24-11-50n.txt';
ABP_filename = 's05237-2926-10-24-11-50_ABP.txt';

estIDs = [1,2,3,4,5,6,7,10,2007];
estNames = ["MAP","Windkessel","SA","SA Warner","Liljestrand","Herd","SA Wesseling","RC Decay","Parlikar"];

figure(10)
set(gcf, 'Position', [100, 100, 1300, 700])
for i = 1 : length(estIDs)
    estID = estIDs(i);
    [feature_times,features_estimated,valsATtco,~] = estimateANDcalibrate(num_filename,ABP_filename,estID);
    
    % C2 calibrate
    CO_estimated = features_estimated(:,1); % C2 calibrated cardiac output [Lpm]
    
    % RMSNE
    TCO_times = valsATtco(:,1); % times of experimental TCO [hr]
    TCO_values = valsATtco(:,2); % values of experimental TCO [Lpm]
    rmsne = RMSNE(TCO_times, TCO_values, feature_times, CO_estimated); % calculate RMSNE
    all_rmsne(i) = rmsne; % each column (=estID) is RMSNE of different estimator
    
    % plot
    subplot(5,2,i)
    plot(feature_times,CO_estimated) % CO from MAP estimator [Lpm]
    hold on;
    stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
    ylim([0 15]);
    title(estNames(i), 'FontSize',12);
end

% Robust Parlikar
[feature_times,features,valsATtco] = RobustParlikar(num_filename,ABP_filename);
% RMSNE
CO_estimated = features(:,1); % cardiac output [Lpm]
TCO_times = valsATtco(:,1); % times of experimental TCO [hr]
TCO_values = valsATtco(:,2); % values of experimental TCO [Lpm]
rmsne_RP = RMSNE(TCO_times, TCO_values, feature_times, CO_estimated); % calculate RMSNE

subplot(5,2,10)
plot(feature_times,CO_estimated) % CO from MAP estimator [Lpm]
hold on;
stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
ylim([0 15]);
xlabel('Time (hr)', 'FontSize',12)
title('Robust Parlikar', 'FontSize',12);

% fix axes labels
subplot(5,2,1); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,3); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,5); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,7); ylabel('CO (Lpm)', 'FontSize',12);
subplot(5,2,9); ylabel('CO (Lpm)', 'FontSize',12); xlabel('Time (hr)', 'FontSize',12);
sgtitle('Patient 5237: Comparison of Different CO Estimation Algorithms')

figure(11)
bar([all_rmsne(1), all_rmsne(2), all_rmsne(3), all_rmsne(4), all_rmsne(5), all_rmsne(6), all_rmsne(7), all_rmsne(8), all_rmsne(9), rmsne_RP])
ylabel('RMSNE', 'FontSize',12);
xticklabels({'MAP','Windkessel','SA','SA Warner','Liljestrand','Herd','SA Wesseling', 'RC Decay', 'Parlikar', 'Robust Parlikar'})
xtickangle(45)
title('Patient 5237: RMSNE Different CO Estimation Algorithms', 'FontSize',12);

%% I: (Q5) Patient 20: Robust Parlikar Estimation of CO, PP, MAP, HR
clear variables; clc;

num_filename = 's00020-2567-03-30-17-47n.txt';
ABP_filename = 's00020-2567-03-30-17-47_ABP.txt';

[feature_times,features,valsATtco] = RobustParlikar(num_filename,ABP_filename);

% get features
CO = features(:,1); % cardiac output [Lpm]
PP = features(:,2); % pulse pressure [mmHg]
MAP = features(:,3); % mean arterial pressure [mmHg]
HR = features(:,4); % heart rate [beats/min]
% TPR = features(:,6); % total peripheral resistance [L/min]

% get values at experimental TCO
TCO_times = valsATtco(:,1); % [hr]
TCO_values = valsATtco(:,2); % [Lpm]
PP_values = valsATtco(:,3); % [mmHg]
MAP_values = valsATtco(:,4); % [mmHg]
HR_values = valsATtco(:,5); % [beats/min]
% TPR_values = valsATtco(:,7); % [L/min]

% plot CO
figure(12)
subplot(4,1,1)
plot(feature_times,CO) % C2 calibrated CO [Lpm]
ylim([1 10])
hold on;
stem(TCO_times,TCO_values,'filled') % TCO values [Lpm]
ylabel('CO', 'FontSize',12);

% plot PP
subplot(4,1,2)
plot(feature_times,PP) % pulse pressure [mmHg]
hold on;
stem(TCO_times,PP_values,'filled') % PP values [mmHg]
ylabel('PP', 'FontSize',12);

% plot MAP
subplot(4,1,3)
plot(feature_times,MAP) % mean arterial pressure [mmHg]
hold on;
stem(TCO_times,MAP_values,'filled') % MAP values [mmHg]
ylabel('MAP', 'FontSize',12);

% plot HR
subplot(4,1,4)
plot(feature_times,HR) % heart rate [beats/min]
hold on;
stem(TCO_times,HR_values,'filled') % HR values [beats/min]
ylabel('HR', 'FontSize',12);
xlabel('time (hours)', 'FontSize',12);
sgtitle('Patient 20: Robust Parlikar Estimation')

%% J: (Q6) TPR estimation
clear variables; clc;

[feature_times_20,features_20,valsATtco_20] = RobustParlikar('s00020-2567-03-30-17-47n.txt','s00020-2567-03-30-17-47_ABP.txt');
[feature_times_214,features_214,valsATtco_214] = RobustParlikar('s00214-3084-11-28-16-23n.txt','s00214-3084-11-28-16-23_ABP.txt');
[feature_times_5237,features_5237,valsATtco_5237] = RobustParlikar('s05237-2926-10-24-11-50n.txt','s05237-2926-10-24-11-50_ABP.txt');

% get TPR feature
TPR_20 = features_20(:,6);
TPR_214 = features_214(:,6);
TPR_5237 = features_5237(:,6);

% plot TPR
figure(13)
subplot(3,1,1)
plot(feature_times_20,TPR_20) % patient 20
ylim([0 4]);
title('Patient 20', 'FontSize',12);

subplot(3,1,2)
plot(feature_times_214,TPR_214) % patient 214
ylim([0 4]);
ylabel('TPR [mmHg/(mL/sec)]', 'FontSize',12);
title('Patient 214', 'FontSize',12);

subplot(3,1,3)
plot(feature_times_5237,TPR_5237) % patient 5237
ylim([0 4]);
xlabel('time (hours)', 'FontSize',12);
title('Patient 5237', 'FontSize',12);
sgtitle('Total Peripheral Resistance Estimations')

%% K: Functions
function [ABP,onsets,EoS_m1,EoS_m2] = abpANDfeatures(ABP_filename, time2start, n)
% ABPANDFEATURES Generates ABP waveform, onsets, and end of systole values
%  In:   ABP_FILENAME   File containing ABP data (eg. 's00020-2567-03-30-17-47_ABP.txt')
%        TIME2START     Time to start ABP wave (in hours)
%        N              Number of ABP pulses 
%
%  Out:  ABP            ABP waveform
%        ONSETS
%               Col 1:  Sample times for onset
%               Col 2:  ABP values for onset
%        EOS_M1
%               Col 1:  Sample times for end of systole using method 1
%               Col 2:  ABP values for end of systole using method 1
%        EOS_M2
%               Col 1:  Sample times for end of systole using method 2
%               Col 2:  ABP values for end of systole using method 2
%
%  method 1: 0.3*sqrt(beatperiod)
%  method 2: lowest non-neg slope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Import and Process Data

% load the ABP data
abp_T = readtable(ABP_filename);
abp_T.Properties.VariableNames = {'Time_sec', 'ABP_mmHg'}; % add headers: time [sec], ABP [mmHg]
times_sec = table2array(abp_T(:,1)); % extract time data from table as an array
times_hrs = times_sec ./ 3600; % convert times from seconds to hours
all_ABP = table2array(abp_T(:,2)); % extract ABP data from table as an array

% first sample when time first above time2start
starting_sample = find(times_hrs >= time2start, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Use 'wabp.m' to get onset times (samples) 

% "Obtains onset sample time (in samples) of each beat in the ABP waveform."
% r = wabp(abp)

all_onsetSamples = wabp(all_ABP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Use 'abpfeature.m' to get features 

% "Extracts features from ABP waveform."
% r = abpfeature(abp,OnsetTimes)

all_features = abpfeature(all_ABP,all_onsetSamples); 

% find end of systole times
all_EoS_m1 = all_features(:,9);  % method 1: 0.3*sqrt(beatperiod)
all_EoS_m2 = all_features(:,11); % method 2: lowest non-neg slope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Onsets

% find samples (x values)
first_onset = find(all_onsetSamples >= starting_sample,1); % ordinal number of first onset after time2start
onsets_n = all_onsetSamples(first_onset:first_onset+(n-1)); % sample time of n onsets after time2start
onsets_samples = (onsets_n - starting_sample); % set first onset to zero

% find ABP values (y values)
onsets_ABP = all_ABP(onsets_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ABP Waveform (num2plot pulses)

ABP = all_ABP(starting_sample:onsets_n(n)); % [mmHg] n ABP pulses after time2start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           End of Systole (method 1)

% find samples (x values)
first_EoS_m1 = find(all_EoS_m1 >= starting_sample,1); % ordinal number of first EoS after time2start
EoS_m1_n = all_EoS_m1(first_EoS_m1:first_EoS_m1+(n-1)); % sample time of n EoS after time2start
EoS_m1_samples = (EoS_m1_n - starting_sample); % set first onset to zero

% find ABP values (y values)
EoS_m1_ABP = all_ABP(EoS_m1_n); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           End of Systole (method 2)

% find samples (x values)
first_EoS_m2 = find(all_EoS_m2 >= starting_sample,1); % ordinal number of first EoS after time2start
EoS_m2_n = all_EoS_m2(first_EoS_m2:first_EoS_m2+(n-1)); % sample time of n EoS after time2start
EoS_m2_samples = (EoS_m2_n - starting_sample); % set first onset to zero

% find ABP values (y values)
EoS_m2_ABP = all_ABP(EoS_m2_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Remove end of systole value if sample time after n'th onset

% end of systole method 1
if EoS_m1_samples(n) > onsets_samples(n)
    EoS_m1_samples = EoS_m1_samples(1:end-1);
    EoS_m1_ABP = EoS_m1_ABP(1:end-1);
end

% end of systole method 2
if EoS_m2_samples(n) > onsets_samples(n)
    EoS_m2_samples = EoS_m2_samples(1:end-1);
    EoS_m2_ABP = EoS_m2_ABP(1:end-1);
end

% final features
onsets = [onsets_samples, onsets_ABP];
EoS_m1 = [EoS_m1_samples, EoS_m1_ABP];
EoS_m2 = [EoS_m2_samples, EoS_m2_ABP];

end

function [C2_vals,k] = C2calibration(TCO_times, TCO_vals, uncal_times, uncal_vals)
% C2CALIBRATION Computes C2 calibration of uncalibrated features
%  In:   TCO_TIMES   	<1xn> times of TCO measurments
%        TCO_values     <1xn> measured TCO values
%        UNCAL_TIMES    <mx1> times of uncalibrated feature measurements
%        UNCAL_VALS     <mx1> uncalibrated feature values
%
%  Out:  C2_VALS        <mx1> C2 calibrated features
%        K              <1x1> calibration factor
%
%  algorithm: first TCO / estimated feature at time of first TCO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

first_TCO_val = TCO_vals(1); % first TCO measurement
first_TCO_time = TCO_times(1); % time of first TCO measurement
first_TCO_idx = find(uncal_times >= first_TCO_time,1);
feature_at_tco = uncal_vals(first_TCO_idx); % estimated feature at time of first TCO
k = first_TCO_val / feature_at_tco; % calibration factor
C2_vals = uncal_vals * k; % C2 calibrated values
end

function [feature_times,features,valsATtco,C] = estimateANDcalibrate(num_filename,ABP_filename,estID)
% ESTIMATEANDCALIBRATE Generates C2 calibrated estimates of CO, PP, MAP, HR
% for the first 12 hours of data
%  NOTE: manually remove row of units (second row) in '*n.txt' file before importing
%  In:   NUM_FILENAME   File containing numerical data ('*n.txt'    -> 's00020-2567-03-30-17-47n.txt')
%        ABP_FILENAME   File containing ABP data       ('*_ABP.txt' -> 's00020-2567-03-30-17-47_ABP.txt')
%        ESTID          CO estimator to use
%                         1:  est01_MAP      - Mean pressure
%                         2:  est02_WK       - Windkessel 1st order LTI RC circuit model
%                         3:  est03_SA       - Systolic area distributed model
%                         4:  est04_SAwarner - Warner systolic area with time correction
%                         5:  est05_Lilj     - Liljestrand PP/(Psys+Pdias) estimator
%                         6:  est06_Herd     - Herd estimator
%                         7:  est07_SAwessCI - Wesseling systolic area with impedance correction
%                         10: est10_RCdecay  - RC exponential decay fit
%                         2007: estimate_parlikar - Parlikar estimator
%
%  Out:  FEATURE_TIMES   <mx1> times of feature measurements
%        FEATURES
%               Col 1:   <mx1> C2 calibrated cardiac output [Lpm]
%               Col 2:   <mx1> pulse pressure [mmHg]
%               Col 3:   <mx1> mean arterial pressure [mmHg]
%               Col 4:   <mx1> heart rate [beats/min]
%               Col 5:   <mx1> total peripheral resistance
%        VALSATTCO
%               Col 1:   <1xn> times of TCO measurments [hr]
%               Col 2:   <1xn> measured TCO values [Lpm]
%               Col 3:   <1xn> PP values when TCO measured [mmHg]
%               Col 4:   <1xn> MAP values when TCO measured [mmHg]
%               Col 5:   <1xn> HR values when TCO measured [beats/min]
%               Col 6:   <1xn> TPR values when TCO measured 
%        C               <1x1> calibration factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the numerical data
num_T = readtable(num_filename);

% find column index for times and TCO
Times_idx = find(string(num_T.Properties.VariableNames) == "ElapsedTime"); % time [sec]
TCO_idx = find(string(num_T.Properties.VariableNames) == "CO"); % thermodilution cardiac output [Lpm]

% extract times and TCO
all_times = table2array(num_T(:,Times_idx)) / 60; % time [min]
all_TCO = table2array(num_T(:,TCO_idx)); % thermodilution cardiac output [Lpm]

% load the ABP data
abp_T = readtable(ABP_filename);
all_ABP = table2array(abp_T(:,2)); % extract ABP data from table as an array

% get data from functions
all_onsetSamples = wabp(all_ABP); % get onset times [samples]
all_features = abpfeature(all_ABP,all_onsetSamples); % get all features
[BeatQ, ~] = jSQI(all_features, all_onsetSamples, all_ABP); % BEATQ: logical of each beat

% estimate CO (uncalibrated)
if estID ~= 2007
    [co, to, ~, fea] = estimateCO_v3(all_onsetSamples,all_features,BeatQ,estID,0);
else
    [co, to, ~, fea] = estimate_parlikar(all_ABP,0);
    delta_P = fea(:,13); % beat-to-beat pressure change at onset times
end
% co - uncalibrated CO [Lpm]
% to - time [min]
% fea - feature matrix
to = to/60; % [hr]

% get features
PP = fea(:,5); % pulse pressure [mmHg]
MAP = fea(:,6); % mean arterial pressure [mmHg]
beat_period = fea(:,7); % beat period [samples]
HR = 60*125./beat_period; % heart rate [beats/min]

% get values when TCO recorded
TCO_times = [];
for i = 1 : length(all_TCO)
    if all_TCO(i) ~= 0
        j = length(TCO_times);
        TCO_times(j+1) = all_times(i) / 60; % [hr]
        TCO_vals(j+1) = all_TCO(i); % [Lpm]
        PP_vals(j+1) = PP(i); % [mmHg]
        MAP_vals(j+1) = MAP(i); % [mmHg]
        HR_vals(j+1) = HR(i); % [beats/min]
    end
end

% C2 calibration
[co_C2cal,k] = C2calibration(TCO_times, TCO_vals, to, co); % cardiac output


% calculate TPR
if estID == 2007
    % MAP / ( CO - k(dPP/beatperiod))
    TPR = MAP./(co_C2cal - k.*(delta_P./beat_period)); % total peripheral resistance

    % get value when TCO recorded
    TPR_vals = [];
    for i = 1 : length(all_TCO)
        if all_TCO(i) ~= 0
            j = length(TPR_vals);
            TPR_vals(j+1) = TPR(i);
        end
    end
end

% get only first 12 hours
TPR_12hrs = [];
for i = 1 : length(to)
    if to(i) <= 12
        to_12hrs(i) = to(i);
        co_12hrs(i) = co_C2cal(i);
        PP_12hrs(i) = PP(i);
        MAP_12hrs(i) = MAP(i);
        HR_12hrs(i) = HR(i);
        if estID == 2007
            TPR_12hrs(i) = TPR(i);
        end
    end
end
TPRv_12hrs = [];
for i = 1 : length(TCO_times)
    if TCO_times(i) <= 12
        TCOt_12hrs(i) = TCO_times(i);
        TCOv_12hrs(i) = TCO_vals(i);
        PPv_12hrs(i) = PP_vals(i);
        MAPv_12hrs(i) = MAP_vals(i);
        HRv_12hrs(i) = HR_vals(i);
        if estID == 2007
            TPRv_12hrs(i) = TPR_vals(i);
        end
    end    
end

% output
feature_times = to_12hrs;
features = [co_12hrs;PP_12hrs;MAP_12hrs;HR_12hrs;TPR_12hrs].';
valsATtco = [TCOt_12hrs;TCOv_12hrs;PPv_12hrs;MAPv_12hrs;HRv_12hrs;TPRv_12hrs].';
C = k; %calibration factor, or Cn in Parkilar
end

function rmsne = RMSNE(TCO_times, TCO_values, CO_times, CO_values)
% RMSNE Calculate Root Mean Squared Normalized Error
%   In:     TCO_times   <1xn> times of TCO measurments [hr]
%           TCO_values  <1xn> measured TCO values [Lpm]
%           CO_times    <mx1> times of feature measurements
%           CO_values   <mx1> C2 calibrated cardiac output [Lpm]
%
%   Out:    RMSNE        <1x1> rmsne value 

% find estimated CO values at times of TCO measurement
for i = 1 : length(TCO_times)
    idx = find(CO_times >= TCO_times(i), 1);
    CO_at_TCO(i,1) = CO_values(idx);
end

% calculate root mean squared normalized error
rmsne = sqrt(mean(((TCO_values - CO_at_TCO)./TCO_values).^2));  
end

function [co, to, told, fea] = estimate_parlikar(abp,filt_order)
% ESTIMATE_PARLIKAR estimate uncalibrated CO using a mothod discussed in
% Parlikar's paper (2007)
%   In:     ABP
%           FILT_ORDER
%   Out:    CO         <kx1>  --- estimated CO (uncalibrated)
%           TO         <kx1>  --- time [minutes] (not evenly sampled!)
%           TOLD       <mx1>  --- time [minutes], not sqi filtered
%         	FEA        <kx13> --- feature matrix

% prepare inputs for estimation
T_Onset = wabp(abp); % sample time of onsets
ABP_Onset = abp(T_Onset); % ABP at onsets

feat = abpfeature(abp, T_Onset);
Pdias = feat(:,4); %Diastolic BP
MAP = feat(:,6); % Mean Pressure (P bar)
Period = feat(:,7); % Period, samples
Period_min = Period/(125*60); % Period, minutes
delta_P = diff(ABP_Onset); % beat-to-beat pressure change at onset times
feat(:,13) = delta_P; % add beat-to-beat pressure change to the feature matrix
[beatq, ~] = jSQI(feat, T_Onset, abp);

% hyperparameter
alpha = 2; % constant, 2 in paper (Parlikar 2007)
L_window = 101; % odd number, slicing window length

%  estimate of 1/τn for every beat, use least-squire-error method
%  consider re_Tau as slope in linear regression with zero intercept
X = MAP;
Y = (alpha*( MAP - Pdias)- delta_P)./Period_min;
re_Tau = zeros(size(Period)); % 1/Tau for each pulse
for i = 1 : length(Period)
    if (i < (L_window+1)/2) || (i > length(Period)-(L_window-1)/2) % leave room for slicing
        continue
    else
        X_s = X(i-(L_window-1)/2:i+(L_window-1)/2);
        Y_s = Y(i-(L_window-1)/2:i+(L_window-1)/2);
        re_Tau(i) = sum(X_s.*Y_s)/sum(X_s.^2);
    end
end

CO_uncalibrated = (delta_P./Period_min + MAP.*re_Tau);

x = CO_uncalibrated;

% codes blow are copied from estimateCO_v3
% synchronize time of of all segments
t = T_Onset(1:end-1)/(125*60); %[min]

% remove all datapoints with bad SQI
sqi1   = beatq(:,1);
ind    = find(sqi1);
x(ind) = [];
t(ind) = [];
t_old  = t;

% apply zero-phase moving avg LPF
if filt_order<2, x_filt = x;
else             x_filt = filtfilt(ones(filt_order,1)/filt_order,1,x);
end

% features
ff = feat;
ff(ind,:) = [];
    
% append segment data to pool
co = x_filt;
to = t;
told = t_old;
fea = ff;
end

function [feature_times,features,valsATtco] = RobustParlikar(num_filename,ABP_filename)
% ROBUSTPARLIKAR implement calibration method described in Parlikar's
% paper, rather than C2 method. It may thus generate a better estimation of
% CO and TPR.
%  In:   NUM_FILENAME   File containing numerical data ('*n.txt'    -> 's00020-2567-03-30-17-47n.txt')
%        ABP_FILENAME   File containing ABP data       ('*_ABP.txt' -> 's00020-2567-03-30-17-47_ABP.txt')
%        
%  Out:  FEATURE_TIMES   <mx1> times of feature measurements [hr]
%        FEATURES        <mx6>
%               Col 1:   <mx1> C2 calibrated cardiac output [Lpm]
%               Col 2:   <mx1> pulse pressure [mmHg]
%               Col 3:   <mx1> mean arterial pressure [mmHg]
%               Col 4:   <mx1> heart rate [beats/min]
%               Col 5:   <mx1> Cn, calibration factor for each estimated CO
%               Col 6:   <mx1> total peripheral resistance [mmHg/(mL/sec)]
%        VALSATTCO       <1x6>
%               xxxxx%Col 1:   <1xn> TCO_times [hr]
%               xxxxx%Col 2:   <1xn> TCO_vals [Lpm]

%               Col 1:   <1xn> times of TCO measurments [hr]
%               Col 2:   <1xn> measured TCO values [Lpm]
%               Col 3:   <1xn> PP values when TCO measured [mmHg]
%               Col 4:   <1xn> MAP values when TCO measured [mmHg]
%               Col 5:   <1xn> HR values when TCO measured [beats/min]
%               Col 6:   <1xn> Cn values when TCO measured
%               Col 7:   <1xn> TPR values when TCO measured [mmHg/(mL/sec)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% <<<<<<<<<<<<<<< hyperparameter >>>>>>>>>>>>>>>>>
alpha = 2; % constant, 2 in paper (Parlikar 2007)
L_window = 21; % odd number, slicing window length
% <<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>

% prepare the recorded TCO data
num_T = readtable(num_filename);
Times_idx = find(string(num_T.Properties.VariableNames) == "ElapsedTime"); % time [sec]
TCO_idx = find(string(num_T.Properties.VariableNames) == "CO"); % thermodilution cardiac output [Lpm]

% extract times and TCO
all_times = table2array(num_T(:,Times_idx)) / 60; % time [min]
all_TCO = table2array(num_T(:,TCO_idx)); % thermodilution cardiac output [Lpm]

% load abp
abp_T = readtable(ABP_filename);
abp = table2array(abp_T(:,2));

% prepare inputs for estimation
T_Onset = wabp(abp); % sample time of onsets
ABP_Onset = abp(T_Onset); % ABP at onsets

feat = abpfeature(abp, T_Onset);
Pdias = feat(:,4); %Diastolic BP
PP = feat(:,5);
MAP = feat(:,6); % Mean Pressure (P bar)
Period = feat(:,7); % Period, samples
Period_min = Period/(125*60); % Period, minutes
HR = 60*125./Period; % heart rate, beats per min
delta_P = diff(ABP_Onset); % beat-to-beat pressure change at onset times
feat(:,13) = delta_P; % add beat-to-beat pressure change to the feature matrix
[beatq, ~] = jSQI(feat, T_Onset, abp);

% get values when TCO recorded
TCO_times = [];
TCO_vals = [];
for i = 1 : length(all_TCO)
    if all_TCO(i) ~= 0
        j = length(TCO_times);
        TCO_times(j+1) = all_times(i) / 60; % [hr]
        TCO_vals(j+1) = all_TCO(i); % [Lpm]
        PP_vals(j+1) = PP(i); % [mmHg]
        MAP_vals(j+1) = MAP(i); % [mmHg]
        HR_vals(j+1) = HR(i); % [beats/min]
    end
end

%  estimate of 1/τn for every beat, use least-squire-error method
%  consider re_Tau as slope in linear regression with zero intercept
X = MAP;
Y = (alpha*( MAP - Pdias)- delta_P)./Period_min;
re_Tau = zeros(size(Period)); % 1/Tau for each pulse
for i = 1 : length(Period)
    if (i < (L_window+1)/2) || (i > length(Period)-(L_window-1)/2) % leave room for slicing
        continue
    else
        X_s = X(i-(L_window-1)/2:i+(L_window-1)/2);
        Y_s = Y(i-(L_window-1)/2:i+(L_window-1)/2);
        re_Tau(i) = sum(X_s.*Y_s)/sum(X_s.^2);
    end
end
Tau = 1./re_Tau;

% compute uncalibrated CO
CO_uncalibrated = (delta_P./Period_min + MAP.*re_Tau);
CO_times = T_Onset(1:end-1)/(125*60*60); %[hr]

% calibration, method mentioned in Parlikar's paper
% calibration factor: Cn
% find the estimated data where there is a corresponding TCO data
index_co = zeros(size(TCO_vals));
for i = 1 : length(TCO_vals)
    index_co(i) = find(CO_times >= TCO_times(i),1);
end
co_corr = CO_uncalibrated(index_co);%not used in computation, just to compare its value with tco
Tau_corr = Tau(index_co);
delta_P_corr = delta_P(index_co);
Period_corr = Period_min(index_co); % min
MAP_corr = MAP(index_co);

% gamma1 and gamma2 computation 
% least square regression
% consider gamma_2 as slope and gamma_1 as intercept
X = MAP_corr;
Y = TCO_vals'./(delta_P_corr./Period_corr + MAP_corr./Tau_corr);
N = length(X);
gamma_2 = ( N*sum(X.*Y) - sum(X)*sum(Y) )/( N*sum(X.^2) - sum(X)^2 );
gamma_1 = (sum(Y) - gamma_2*sum(X)) / N;

% CO calibration using Cn
Cn = gamma_1 + gamma_2.*MAP;
CO_calibrated = Cn .* CO_uncalibrated;

% compute TPR
TPR = MAP ./ (CO_calibrated - Cn.*delta_P./Period_min); % [mmHg/(L/min)]
TPR = TPR .* (60/1000); % [mmHg/(mL/sec)]

% get value when TCO recorded
TPR_vals = [];
for i = 1 : length(all_TCO)
    if all_TCO(i) ~= 0
        j = length(TPR_vals);
        TPR_vals(j+1) = TPR(i);
        Cn_vals(j+1) = Cn(i);
    end
end

% remove bad points
sqi1   = beatq(:,1);
ind    = find(sqi1);
CO_times(ind) = []; 
CO_calibrated(ind) = [];
PP(ind) = [];
MAP(ind) = [];
HR(ind) = [];
Cn(ind) = [];
TPR(ind) = [];

% get only first 12 hours
TPR_12hrs = [];
for i = 1 : length(CO_times)
    if CO_times(i) <= 12
        to_12hrs(i) = CO_times(i);
        co_12hrs(i) = CO_calibrated(i);
        PP_12hrs(i) = PP(i);
        MAP_12hrs(i) = MAP(i);
        HR_12hrs(i) = HR(i);
        Cn_12hrs(i) = Cn(i);
        TPR_12hrs(i) = TPR(i);
    end
end
TPRv_12hrs = [];
for i = 1 : length(TCO_times)
    if TCO_times(i) <= 12
        TCOt_12hrs(i) = TCO_times(i);
        TCOv_12hrs(i) = TCO_vals(i);
        PPv_12hrs(i) = PP_vals(i);
        MAPv_12hrs(i) = MAP_vals(i);
        HRv_12hrs(i) = HR_vals(i);
        Cnv_12hrs(i) = Cn_vals(i);
        TPRv_12hrs(i) = TPR_vals(i);
    end    
end

% output
feature_times = to_12hrs;
features = [co_12hrs;PP_12hrs;MAP_12hrs;HR_12hrs;Cn_12hrs;TPR_12hrs].';
valsATtco = [TCOt_12hrs;TCOv_12hrs;PPv_12hrs;MAPv_12hrs;HRv_12hrs;Cnv_12hrs;TPRv_12hrs].';
end

