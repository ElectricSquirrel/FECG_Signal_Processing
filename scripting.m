% scripting

%% 'walking artifacts' and 'QRS-only 2-lead ECG signals'
% Section1: testing between 
% 'walking artifacts' and 'QRS-only 2-lead ECG signals'
% clear
close all
clc
[t1,exp_data] = rdsamp('mitdb/100', 1, 11000, 201); 
plot(t1,exp_data);
xlim([1 10])
%%
load('test11_45wm.mat')
[tm,signal,Fs,labels]=rdmat('test11_45wm');

% load qrs-only signal, because it's a two lead signal, there are only two
% columns in the signal_qrs data
load('sel100m.mat');
[tm_qrs,signal_qrs,Fs_qrs,labels_qrs]=rdmat('sel100m');


% Plot and compare
subplot(211);plot(tm(1:2000), signal(1:2000,4));
title('ECG walking artifacts singal lead1');xlabel('time domain(s)')
subplot(212);plot(tm_qrs(1:1000), signal_qrs(1:1000,1));
title('QRS only ECG signal lead1')
xlabel('time domain(s)')

%%

% wfdb2mat('macecgdb/test20_45j.dat');
% [tm,signal,Fs,labels]=rdmat('test20_45jm');

%% Fetal Heart Rate Database signals

wfdb2mat('adfecgdb/r10.edf')
[tm,signal,Fs,labels]=rdmat('r10_edfm');


% Column 1: Direct_1
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
suptitle('Sample FHR signals from MIT-BIH arrhythmia database')

% subplot(212);plot(tm_qrs(1:1000), signal_qrs(1:1000,1));




%% Testing for Peyre Code:
% clear
clc
close all
clear

N = 2048*2;
name = 'piece-regular';
f0 = load_signal(name, N);
f0 = rescale(f0,.05,.95);
sigma = 0.05;

f = f0 + randn(size(f0))*sigma;

% plot(1:length(f),f)
clf;
subplot(2,1,1);
plot(f0); axis([1 N 0 1]);
title('Clean signal');
subplot(2,1,2);
plot(f); axis([1 N 0 1]);
title('Noisy signal');


% Hard Thresholding code
% Functionality: If x is smaller than threshold T, then 
% the corresponding signal will be zerod
Theta0 = @(x,T)x .* (abs(x)>T);

% Soft Thresholding code
Theta1 = @(x,T)max(0, 1-T./max(abs(x),1e-9)) .* x;




t = linspace(-3,3,1024)'; T = 1;
clf;
plot( t, [Theta0(t,T), Theta1(t,T)], 'LineWidth', 2 );
axis('equal'); axis('tight');
legend('\Theta^0', '\Theta^1');

options.ti = 0; 
Jmin = 4;
W  = @(f) perform_wavelet_transf(f,Jmin,+1,options);    %figure out what the function is doing 10pm
Wi = @(fw)perform_wavelet_transf(fw,Jmin,-1,options);

x = W(f);
x1 = Theta0(x, 3*sigma);
clf;
subplot(2,1,1);
plot_wavelet(x,Jmin); axis([1 N -1 1]);
title('W(f)');
subplot(2,1,2);
plot_wavelet(Theta0(W(f),T),Jmin); axis([1 N -1 1]);    % apply this 
title('\Theta_0(W(f))');




%%
close all
clc
Theta0 = @(x,T)x .* (abs(x)>T);
t = 1:4000;
T = 0.5;

[swa,swd] = swt(signal(:,2),3,'db1');


swa(1,:) = Theta0(swa(1,:),2)
swa(2,:) = Theta0(swa(2,:),2)
swa(3,:) = Theta0(swa(3,:),2)

swd(1,:) = Theta0(swd(1,:),T)
swd(2,:) = Theta0(swd(2,:),T)
swd(3,:) = Theta0(swd(3,:),T)

subplot(311);plot(t,swa(1,:))
subplot(312);plot(t,swa(2,:))
subplot(313);plot(t,swa(3,:))
% subplot(414);plot(1:4000,swa(4,:))


I_signal = iswt(swa,swd,'db1');
figure;
subplot(211);plot(1:4000,signal(:,2))
title('Original Signal')
subplot(212);plot(1:4000,I_signal)
title('SWT-denoised Signal')

%% Real Data Processing
%
clear
close all
clc
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('wfdb-app-toolbox-0-9-9/mcode')
addpath('toolbox_wavelet_meshes')
Theta0 = @(x,T)x .* (abs(x)>T);

undisturbed = csvread('undisturbed.txt',1,0);
undisturbed(:,1) = undisturbed(:,1)/1000;
t = undisturbed(:,1)
T = 100;

[swa,swd] = swt(undisturbed(:,2)',3,'db1');

swa(1,:) = Theta0(swa(1,:),T)
swa(2,:) = Theta0(swa(2,:),T)
% swa(3,:) = Theta0(swa(3,:),T)

subplot(311);plot(t,swa(1,:))
subplot(312);plot(t,swa(2,:))
subplot(313);plot(t,swa(3,:))
suptitle('SWA (stationary wavelet approximation plot)')

swd(1,:) = Theta0(swd(1,:),T)
swd(2,:) = Theta0(swd(2,:),T)
% swd(3,:) = Theta0(swd(3,:),T)
figure
subplot(311);plot(t,swd(1,:))
subplot(312);plot(t,swd(2,:))
subplot(313);plot(t,swd(3,:))
suptitle('SWD (stationary wavelet Coefficient plot)')

I_signal = iswt(swa,swd,'db1');
figure;
subplot(211);plot(undisturbed(:,1),undisturbed(:,2))
title('Original Signal')
subplot(212);plot(undisturbed(:,1),I_signal)
title('SWT-denoised Signal');
xlabel('time: seconds')

%% wireSwingTxt
%
clear
close all
clc
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('wfdb-app-toolbox-0-9-9/mcode')
addpath('toolbox_wavelet_meshes')
Theta0 = @(x,T)x .* (abs(x)>T);

undisturbed = csvread('flexing30secs.txt',1,0);
undisturbed(:,1) = undisturbed(:,1)/1000;
t = undisturbed(:,1)
T = 100;


[swa,swd] = swt(undisturbed(:,2)',3,'db1');

% swa(1,:) = Theta0(swa(1,:),T)
swa(2,:) = Theta0(swa(2,:),T);
swa(3,:) = Theta0(swa(3,:),T);

subplot(311);plot(t,swa(1,:))
subplot(312);plot(t,swa(2,:))
subplot(313);plot(t,swa(3,:))
suptitle('SWA (stationary wavelet approximation plot)')

% swd(1,:) = Theta0(swd(1,:),T)
swd(2,:) = Theta0(swd(2,:),T)
swd(3,:) = Theta0(swd(3,:),T)
figure
subplot(311);plot(t,swd(1,:))
subplot(312);plot(t,swd(2,:))
subplot(313);plot(t,swd(3,:))
suptitle('SWD (stationary wavelet Coefficient plot)')

I_signal = iswt(swa,swd,'db1');
figure;
subplot(211);plot(undisturbed(:,1),undisturbed(:,2))
title('flexing30secs')
subplot(212);plot(undisturbed(:,1),I_signal)
title('SWT-denoised Signal');
xlabel('time: seconds')