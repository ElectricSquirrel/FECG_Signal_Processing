function sig2 = LPF_butterworth(sig1, Fcut)

% this creates a butter worth filter and cleans out the high frequency noise
% over 20 Hz, including the DC noise

filter_order = 5;

[a,b] = butter(filter_order,Fcut/500);

sig2 = filtfilt(a,b,sig1);


end