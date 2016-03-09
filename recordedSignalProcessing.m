
%% Recorded Signals from Spark fun
% Problem is that the sampling rate is not high enough
% sol1: record a new set of data with higher sampling rate
% sol2: built a new set of 

undisturbed = csvread('undisturbed.txt',1,0);
undisturbed(:,1) = undisturbed(:,1)/1000;


flexing30secs = csvread('flexing30secs.txt',1,0);
flexing30secs(:,1) = flexing30secs(:,1)/1000; %setting the timeframe to seconds

wireSwing30s = csvread('wireSwing30s.txt',1,0);
wireSwing30s(:,1) = wireSwing30s(:,1)/1000;


subplot(311);plot(undisturbed(:,1),undisturbed(:,2))
title('undisturbed Signal')

subplot(312);plot(flexing30secs(:,1),flexing30secs(:,2))
title('flexing30secs Signal')

subplot(313);plot(wireSwing30s(:,1),wireSwing30s(:,2))
title('wireSwing30s Signal')

%% FFT of all three signals

undisturbedFFT = fft(undisturbed(:,2) );

flexing30secsFFT = fft(flexing30secs(:,2) );

wireSwing30sFFT = fft(wireSwing30s(:,2) );

figure

subplot(311);plot(1:length(undisturbedFFT),undisturbedFFT);
title('undisturbed Signal');

subplot(312);plot(1:length(flexing30secsFFT),flexing30secsFFT);
title('flexing30secs Signal');

subplot(313);plot(1:length(wireSwing30sFFT),wireSwing30sFFT);
title('wireSwing30s Signal');
suptitle('FFT of all three signals');
%% Figuring out the SWT frequency correspondence
% swt:
% swt(1) is high frequency
% swt(end) is lowest frequency range
close all
clear
clc
Theta0 = @(x,T)x .* (abs(x)>T);

undisturbed = csvread('undisturbed.txt',1,0);
undisturbed(:,1) = undisturbed(:,1)/1000;

T = 1000;
t = undisturbed(:,1);
T2 = T/10;

[swa,swd] = swt(undisturbed(:,2)',4,'db1');

swa(1,:) = Theta0(swa(1,:),T);
% swa(2,:) = Theta0(swa(2,:),T2);
% swa(3,:) = Theta0(swa(3,:),T);

% subplot(311);plot(t,swa(1,:));
% subplot(312);plot(t,swa(2,:));
% subplot(313);plot(t,swa(3,:));
% suptitle('SWA (stationary wavelet approximation plot)')

swd(1,:) = Theta0(swd(1,:),T);
% swd(2,:) = Theta0(swd(2,:),T2);
% swd(3,:) = Theta0(swd(3,:),T);
% figure
% subplot(311);plot(t,swd(1,:))
% subplot(312);plot(t,swd(2,:))
% subplot(313);plot(t,swd(3,:))
% suptitle('SWD (stationary wavelet Coefficient plot)')


I_signal = iswt(swa,swd,'db1');

undisturbedFFT = fft(undisturbed(:,2) );
I_signalFFT = fft(I_signal);

subplot(211);plot(1:length(undisturbedFFT),undisturbedFFT);
title('Original undisturbed Signal FFT');
subplot(212);plot(1:length(I_signalFFT),I_signalFFT);
title('denoised Signal FFT');

figure
subplot(211);plot(undisturbed(:,1),undisturbed(:,2))
title('undisturbed')
subplot(212);plot(undisturbed(:,1),I_signal)
title('SWT-denoised Signal');
xlabel('time: seconds')
suptitle('Time Domain Plot of two signals');

%% Apply FFT to MITBIT Fetal Heart Rate Database signals

wfdb2mat('adfecgdb/r10.edf')
[tm,signal,Fs,labels]=rdmat('r10_edfm');

% Column 1: Direct_1 on fetus head
% Column 2: Abdomen_1
% Column 3: Abdomen_2
% Column 4: Abdomen_3
% Column 5: Abdomen_4
% Column 6: EDF Annotations


subplot(511);plot(tm(1:4000), signal(1:4000,1));
subplot(512);plot(tm(1:4000), signal(1:4000,2));
subplot(513);plot(tm(1:4000), signal(1:4000,3));
subplot(514);plot(tm(1:4000), signal(1:4000,4));
subplot(515);plot(tm(1:4000), signal(1:4000,5));
suptitle('Sample FHR signal from MIT-BIH arrhythmia database')

% subplot(212);plot(tm_qrs(1:1000), signal_qrs(1:1000,1));

%%
% Applying Butterworth filter to clean out the baseline wander
start_point = 8001;

[a,b] = butter(5,3/500,'high');
% sos = zp2sos(z,p,k);
% fvtool(sos,'Analysis','freq')
sig1 = signal(:,5);
sig2 = filtfilt(a,b,sig1);
subplot(211);plot(1:length(sig1), sig1);
xlim([start_point start_point+4000])
subplot(212);plot(1:length(sig2), sig2);
xlim([start_point start_point+4000])
%% 
% we are applying swt high frequency hard thresholding on 3 and 5

Theta0 = @(x,T)x .* (abs(x)>T);
FECG = sig2;
% FECG = sig2;
t = 1:length(FECG);

T = 1000;
% T2 = T/10;

[swa,swd] = swt(FECG',5,'db1');

swa(1,:) = Theta0(swa(1,:),T);
swa(2,:) = Theta0(swa(2,:),T);
swa(3,:) = Theta0(swa(3,:),T);
swa(4,:) = Theta0(swa(4,:),T);

swd(1,:) = Theta0(swd(1,:),T);
swd(2,:) = Theta0(swd(2,:),T);
swd(3,:) = Theta0(swd(3,:),T);
swd(4,:) = Theta0(swd(4,:),T);

FECG_Denoised = iswt(swa,swd,'db1');

FECG_FFT = fft(FECG);
FECG_Denoised_FFT = fft(FECG_Denoised);

subplot(211);plot(1:length(FECG_FFT),FECG_FFT);
title('Original FECG Signal FFT');
subplot(212);plot(1:length(FECG_Denoised_FFT),FECG_Denoised_FFT);
title('denoised FECG FFT');

figure
subplot(211);plot(t(start_point:start_point+4000),FECG(start_point:start_point+4000))
title('Original FECG Signal')
subplot(212);plot(t(start_point:start_point+4000),FECG_Denoised(start_point:start_point+4000))
title('SWT-denoised Signal');
% xlabel('time: seconds')
suptitle('Time Domain Plot of two signals');

%%
fft_check(150,FECG,'FECG baseline fft');