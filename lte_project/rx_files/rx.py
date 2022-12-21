# Imports
from numpy import arange, exp, pi, floor, complex64, array

# Dumb setup stuff because the API is not set up well and I don't know how to fix it.
# You probably have to create uhd_params and put the path in
import sys
from uhd_params import path_to_uhd_module as path

IP_radio = "192.168.0.165"

sys.path.insert(0, path)
import uhd


usrp = uhd.usrp.MultiUSRP(f'IPAddress={IP_radio}')
# Note: num_samps = rate * 0.001 seconds for 1 ms of output
samples = usrp.recv_num_samps(num_samps=61440, freq=3700e6, rate=61.44e6, channels=[0], gain=10) # units: N, Hz, Hz, list of channel IDs, dB

with open("rx_data", "wb") as out_file:
    samples.tofile(out_file)




