# # Imports
# from numpy import arange, exp, pi, floor, complex64, array

# Dumb setup stuff because the API is not set up well and I don't know how to fix it.
# You probably have to create uhd_params and put the path in
import sys
from uhd_params import path_to_uhd_module as path

IP_radio = "10.2.118.69"

sys.path.insert(0, path)
# import uhd


# usrp = uhd.usrp.MultiUSRP(f'IPAddress={IP_radio}, recv_frame_size=1472')
# # Note: num_samps = rate * 0.001 seconds for 1 ms of output
# samples = usrp.recv_num_samps(num_samps=61440, freq=3700e6, rate=61.44e6, channels=[0], gain=10) # units: N, Hz, Hz, list of channel IDs, dB

# with open("rx_data", "wb") as out_file:
#     samples.tofile(out_file)

import uhd
import numpy as np

usrp = uhd.usrp.MultiUSRP()

num_samps = 61440 # number of samples received
center_freq = 3700e6 # Hz
sample_rate = 61.44e6 # Hz
gain = 10 # dB

usrp.set_rx_rate(sample_rate, 0)
usrp.set_rx_freq(uhd.libpyuhd.types.tune_request(center_freq), 0)
usrp.set_rx_gain(gain, 0)

# Set up the stream and receive buffer
st_args = uhd.usrp.StreamArgs("fc32", "sc16")
st_args.channels = [0]
# st_args.spp = "2000";   # Also tried 200 and 365 here
metadata = uhd.types.RXMetadata()
streamer = usrp.get_rx_stream(st_args)
print(streamer.get_max_num_samps())
recv_buffer = np.zeros((1, 1000), dtype=np.complex64)

# Start Stream
stream_cmd = uhd.types.StreamCMD(uhd.types.StreamMode.start_cont)
stream_cmd.stream_now = True
streamer.issue_stream_cmd(stream_cmd)

# Receive Samples
samples = np.zeros(num_samps, dtype=np.complex64)
for i in range(num_samps//1000):
    streamer.recv(recv_buffer, metadata)
    samples[i*1000:(i+1)*1000] = recv_buffer[0]

# Stop Stream
stream_cmd = uhd.types.StreamCMD(uhd.types.StreamMode.stop_cont)
streamer.issue_stream_cmd(stream_cmd)

print(len(samples))
print(samples[0:10])




