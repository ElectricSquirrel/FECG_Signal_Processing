function sig2 = LPF_butterworth(sig1)
% this creates a butter worth filter and cleans out the low frequency noise
% below 3 Hz, including the DC noise

filter_order = 5;

[a,b] = butter(filter_order,3/500,'high');

sig2 = filtfilt(a,b,sig1);


end