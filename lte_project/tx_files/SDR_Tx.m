%% Initialize
fc = 3700e6;
ts = 1/(39e6);
load soundingSignal40-39.mat
load ../../itc_paper/matlab/soundingSignal-BW10-Fs20.mat
x_tilde_once = s.'/10;
x_tilde = repmat(x_tilde_once, [1e5 1]);

if ~(exist('tx', 'var'))
    tx = comm.SDRuTransmitter(              ...
                  'Platform','B210',        ...
                  'SerialNum','32339F7',    ...
                  'CenterFrequency',fc,     ...
                  'MasterClockRate',1/ts,   ...
                  'InterpolationFactor',1);
end
tx.Gain = 40;

%% Transmit until stopped
info(tx)
disp('Transmitting')
while 1
    tx(x_tilde);
end 

%% Release the tx (must be run after stopping)
release(tx);

