# Imports
from numpy import arange, exp, pi, floor, complex64, array

# Dumb setup stuff because the API is not set up well and I don't know how to fix it.
# You probably have to create uhd_params and put the path in
import sys
from uhd_params import path_to_uhd_module as path

IP_radio = "192.168.0.165"

sys.path.insert(0, path)
import uhd


usrp = uhd.usrp.MultiUSRP('IPAddress=192.168.0.165')
samples = usrp.recv_num_samps(num_samps=12000000, freq=3700e6, rate=48e6, channels=[0], gain=10) # units: N, Hz, Hz, list of channel IDs, dB
# samples = usrp.recv_num_samps(fc=3700e6, duration=1, sample_rate=24e6, channels=[0], gain=1)

print(samples)




