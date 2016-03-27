
clear
clc
close all

% wfdb2mat('tpehgdb/tpehg1132');
% [tm,signal,Fs,labels]=rdsamp('tpehgdb/tpehg1021');
%% tpehgdb
% Each database system
% Each record is composed of three channels, recorded from 4 electrodes:
% the first electrode (E1) was placed 3.5 cm to the left and 3.5 cm above the navel;
% the second electrode (E2) was placed 3.5 cm to the right and 3.5 cm above the navel;
% the third electrode (E3) was placed 3.5 cm to the right and 3.5 cm below the navel;
% the fourth electrode (E4) was placed 3.5 cm to the left and 3.5 cm below the navel.

% Electrodes Placement
%
% 1. first channel, unfiltered (S1);
% 2. first channel, (S1) filtered using a 4-pole band-pass Butterworth filter from 0.08Hz to 4Hz;
% 3. first channel, (S1) filtered using a 4-pole band-pass Butterworth filter from 0.3Hz to 3Hz;
% 4. first channel, (S1) filtered using a 4-pole band-pass Butterworth filter from 0.3Hz to 4Hz;
% 5. second channel, unfiltered (S2);
% 6. second channel, (S2) filtered using a 4-pole band-pass Butterworth filter from 0.08Hz to 4Hz;
% 7. second channel, (S2) filtered using a 4-pole band-pass Butterworth filter from 0.3Hz to 3Hz;
% 8. second channel, (S2) filtered using a 4-pole band-pass Butterworth filter from 0.3Hz to 4Hz;
% 9. third channel, unfiltered (S3);
% 10. third channel, (S3) filtered using a 4-pole band-pass Butterworth filter from 0.08Hz to 4Hz;
% 11. third channel, (S3) filtered using a 4-pole band-pass Butterworth filter from 0.3Hz to 3Hz;
% 12. third channel, (S3) filtered using a 4-pole band-pass Butterworth filter from 0.3Hz to 4Hz.

[t1,exp_data] = rdsamp('tpehgdb/tpehgdb/tpehg1022',[],35200,200); 
t2 = t1./60; % t1 is in unit seconds, t2 is in unit minute

%%
close all
cp_select = 3;
plot(t2,exp_data(:,[cp_select cp_select+4 cp_select+8]))
title('Term-Preterm EHG Database')
% 