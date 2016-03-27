%%
% This scripting file is for Fetal Heart Rate Extraction Simulation
% 1) included a 2-lead ECG signal from from MIT-BIH databse
% 2) manipulate signal at first channel:
%    High-pass filtered by 1.5Hz
%    downsample by 2 (making it looking like fetal signal, 
%    almost twice as fast as maternal ECG). Add both together.
% 3) PCA: use SWT to clean out the high frequency noise subtract
%    the maternal signal (channel 2 signal) from the combined 
%    signal
% 4) Extract the frequency information of FECG
clear
close all
clc
[t1,exp_data,Fs] = rdsamp('qtdb/sel100'); 
%% Step2

MECG = exp_data(:,2);
MECG = HPF_butterworth(MECG);
FECG = downsample(exp_data(:,1),2);
FECG = HPF_butterworth(FECG);
% Fcut = 400;
% FECG2 = LPF_butterworth(FECG,Fcut);
ComboECG = FECG +MECG(1:length(FECG));

miu = 0.1;  % scale factor
FECG = FECG * miu;
t2 = t1(1:length(ComboECG));

plot(t2,ComboECG)
% plot(t2,FECG2,'-b',t2,FECG,'-r');
% xlim([1 10])

% plot();%,t2,MECG(1:length(FECG)));

xlim([1 5])


%% Applying SWT to the Combined ECG signal
start_point = 1;
Theta0 = @(x,T)x .* (abs(x)>T);

% the line below should only run one time because it is to cope with the
% SWT length requirement
ComboECG = [ComboECG; zeros(12,1)];
t = t1(1:length(ComboECG));

T = 1000;
[swa,swd] = swt(ComboECG',5,'db1');

swa(1,:) = Theta0(swa(1,:),T);
swa(2,:) = Theta0(swa(2,:),T);
% swa(3,:) = Theta0(swa(3,:),T);
% swa(4,:) = Theta0(swa(4,:),T);

swd(1,:) = Theta0(swd(1,:),T);
swd(2,:) = Theta0(swd(2,:),T);
% swd(3,:) = Theta0(swd(3,:),T);
% swd(4,:) = Theta0(swd(4,:),T);

ComboECG_Denoised = iswt(swa,swd,'db1');



figure
subplot(211);plot(t(start_point:start_point+1000),ComboECG(start_point:start_point+1000))
title('Original ComboECG Signal')
subplot(212);plot(t(start_point:start_point+1000),ComboECG_Denoised(start_point:start_point+1000))
title('SWT-denoised Signal');
xlabel('time: seconds')
suptitle('Time Domain Plot of two signals');


%%



MECG = MECG(1:112512);
t = t1(1:length(MECG));

T = 1000;
[swa_M,swd_M] = swt(MECG',5,'db1');

swa_M(1,:) = Theta0(swa_M(1,:),T);
swa_M(2,:) = Theta0(swa_M(2,:),T);
% swa_M(3,:) = Theta0(swa_M(3,:),T);
% swa_M(4,:) = Theta0(swa_M(4,:),T);

swd_M(1,:) = Theta0(swd_M(1,:),T);
swd_M(2,:) = Theta0(swd_M(2,:),T);
% swd_M(3,:) = Theta0(swd_M(3,:),T);
% swd_M(4,:) = Theta0(swd_M(4,:),T);

MECG_Denoised = iswt(swa_M,swd_M,'db1');


% figure
subplot(211);plot(t(start_point:start_point+1000),MECG(start_point:start_point+1000))
title('Original MECG Signal')
subplot(212);plot(t(start_point:start_point+1000),MECG_Denoised(start_point:start_point+1000))
title('SWT-denoised Signal');
xlabel('time: seconds')
suptitle('Time Domain Plot of two signals');


%% Step3: Applying PCA to SWT'ed ComboECG and MECG
F_PCAed_ECG = (ComboECG_Denoised-MECG_Denoised);

subplot(311);plot(t(start_point:start_point+1000),ComboECG_Denoised(start_point:start_point+1000))
ylim([-0.3 0.9]);
title('Denoised ComboECG Signal')

subplot(312);plot(t(start_point:start_point+1000),MECG_Denoised(start_point:start_point+1000))
title('Denoised MECG Signal');
ylim([-0.3 0.9]);

subplot(313);plot(t(start_point:start_point+1000),F_PCAed_ECG(start_point:start_point+1000))
title('FECG_PACed Signal');
ylim([-0.3 0.9]);

xlabel('time: seconds')
suptitle('Time Domain Plot of two signals');

%% Step4: Extract the frequency information
close all
plot(t,F_PCAed_ECG)%,'ro');
xlim([0 5])
% close all

Fs = 1/0.004;
%%
% qrs_i_raw is the index of peak heart rate in "t"
% contains the heart rate time lapse information. We are going to
% use this to extract the FHR RR interval
[qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(F_PCAed_ECG,Fs,1);

%%
%
clc
counter = 5;
temp = zeros(1,counter);
qrs_i_rand_start = randi(size(qrs_i_raw)-counter)
for i = 1:counter
    temp(i) = t(qrs_i_raw(qrs_i_rand_start+i+1))- t(qrs_i_raw(qrs_i_rand_start+i));
%     temp(i) = t(qrs_i_raw(i+1))- t(qrs_i_raw(i));
end
% RR_Tavg is the average RR interval over a "counter" amount of samples
% The unit is "second". It should be below 1 second. Normally around 0.5s
RR_Tavg = sum(temp)/counter

% Fetal Heart Rate per minute, physiological FHR should be around 120 to
% 150 beats per minute
FHR_per_min = 60/RR_Tavg

%%


