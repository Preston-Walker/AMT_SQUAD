% run this one
plutoradiosetup();
configurePlutoRadio('AD9364');

tx = sdrtx('Pluto');
tx.RadioID = 'usb:0';
tx.CenterFrequency = 3.7e9;
tx.BasebandSampleRate = 48e6;
tx.Gain = 0;

load soundingSignal40-48.mat;

% run this every time 

tx.transmitRepeat(s.');