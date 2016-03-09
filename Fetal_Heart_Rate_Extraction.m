%%
% This scripting file is for Fetal Heart Rate Extraction Simulation
% 1) included a 2-lead ECG signal from from MIT-BIH databse
% 2) manipulate signal at first channel:
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
MECG = LPF_butterworth(MECG);
FECG = downsample(exp_data(:,1),2);
FECG = LPF_butterworth(FECG);
% Fcut = 400;
% FECG2 = HPF_butterworth(FECG,Fcut);
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
%ComboECG = [ComboECG; zeros(12,1)];
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



MECG = MECG(1:112512)
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
title('ComboECG - MECG Signal');
ylim([-0.3 0.9]);

xlabel('time: seconds')
suptitle('Time Domain Plot of two signals');

%% Step4: Extract the frequency information

F_PCAed_ECG_FFT = fft(F_PCAed_ECG);

figure
plot(1:length(F_PCAed_ECG_FFT),F_PCAed_ECG_FFT);
title('denoised F-PCAed-ECG FFT');


