"""
Contains our frame synchronization SM
"""

from numpy import mean, argmax, genfromtxt, array, append, dtype, empty, shape, prod
from enum import Enum
#from buffer import FS_buffer
from params import *
from FS_functions import *
from matplotlib.pyplot import plot, text, grid, show, subplot, figure, subplots

st = Enum(
    "FSSM_st",
    "init coldStart onePeak twoPeaks weird locked oneMissed twoMissed",
)


class FSSM:
    def __init__(self, corr="L0_single"):
        self.state = st.init
        self.corr_data = zeros(256)
        self.current_idx = -1
        self.last_mhat = None
        self.llast_mhat = None
        self.mhats = array([])
        self.weird = None
        self.last_weird = None
        self.correlations = array([])
        self.all_correlated_samples = array([])
        self.window_locked = array([])
        self.lock_pairs = array([[]])
        self.lock_begin = None
        self.lock_end = None
        if corr == "L0_single":
            self.corr = L0_single
        elif corr == "L6_single":
            self.corr = L6_single
        else:
            print("Unsupported correlation function")

    def get_peak(self):
        self.correlations = append(self.correlations, (self.corr(self.corr_data)))
        self.correlations = self.correlations[1:]
        avg_val = mean(self.correlations)
        max_idx = argmax(self.correlations)
        max_val = self.correlations[max_idx]
        #print(argmax(self.correlations)+self.current_idx-buffer_size)
        return avg_val, max_idx, max_val

    def tick(self, newSample):
        self.corr_data = append(self.corr_data, newSample)
        self.corr_data = self.corr_data[1:]
        self.all_correlated_samples = append(self.all_correlated_samples, (self.corr(self.corr_data)))
        self.current_idx += 1

        # State update first
        if self.state == st.init:
            if self.correlations.size < buffer_size-256:
                self.state = st.init
                self.correlations = append(self.correlations, (self.corr(self.corr_data)))
            else:
                self.state = st.coldStart
                self.correlations = append(self.correlations, (self.corr(self.corr_data)))
        elif self.state == st.coldStart:
            avg_val, max_idx, max_val = self.get_peak()
            tmp_idx = self.current_idx - buffer_size + max_idx
            if max_val > 3 * avg_val:                              
                self.state = st.onePeak
                self.last_mhat = tmp_idx
                self.mhats = append(self.mhats, self.last_mhat)
            else:
                self.state = st.coldStart
        elif self.state == st.onePeak:
            avg_val, max_idx, max_val = self.get_peak()
            tmp_idx = self.current_idx - buffer_size + max_idx
            if tmp_idx == self.last_mhat:
                self.state = st.onePeak
            elif abs(tmp_idx - (self.last_mhat + frame_size)) <= 1:
                self.state = st.twoPeaks
                self.llast_mhat = self.last_mhat
                self.last_mhat = tmp_idx
                self.mhats = append(self.mhats, self.last_mhat)
            elif abs(tmp_idx - (self.last_weird + frame_size)) <= 1:
                self.state = st.twoPeaks
                self.llast_mhat = self.last_weird
                self.last_mhat = tmp_idx
                self.mhats = append(self.mhats, self.last_mhat)
                self.mhats = append(self.mhats, self.llast_mhat)
            else:
                self.state = st.weird                                       
                self.weird = self.current_idx - buffer_size + max_idx
        elif self.state == st.twoPeaks:
            avg_val, max_idx, max_val = self.get_peak()
            tmp_idx = self.current_idx - buffer_size + max_idx
            if tmp_idx == self.last_mhat:
                self.state = st.twoPeaks
            elif abs(tmp_idx - (self.last_mhat + frame_size)) <= 1:
                self.state = st.locked
                self.llast_mhat = self.last_mhat
                self.last_weird = None
                self.last_mhat = tmp_idx
                self.mhats = append(self.mhats, self.last_mhat)
                self.lock_begin = self.current_idx
            elif abs(tmp_idx - (self.llast_mhat + frame_size)) <= 1:
                self.state = st.twoPeaks
                self.last_mhat = tmp_idx
                self.mhats = append(self.mhats, self.last_mhat)
            else:
                self.state = st.weird
                self.weird = tmp_idx
        elif self.state == st.weird:
            if self.current_idx < self.weird + frame_size:
                self.state == st.weird
            else:
                avg_val, max_idx, max_val = self.get_peak()
                tmp_idx = self.current_idx - buffer_size + max_idx
                if abs(tmp_idx - (self.last_mhat + frame_size)) <= 1:
                    self.state = st.twoPeaks if not self.llast_mhat else st.locked
                    self.llast_mhat = self.last_mhat
                    self.last_mhat = tmp_idx
                    self.mhats = append(self.mhats, self.last_mhat)
                elif abs(tmp_idx - (self.weird + frame_size)) <= 1:
                    self.state = st.twoPeaks
                    self.llast_mhat = self.weird
                    self.last_mhat = tmp_idx
                    self.mhats = append(self.mhats, self.last_mhat)
                else:
                    self.state = st.onePeak
                    self.last_weird = self.weird
                    self.last_mhat = tmp_idx
                    self.mhats = append(self.mhats, self.last_mhat)
                    self.llast_mhat = None
                self.weird = None
        elif self.state == st.locked:
            avg_val, max_idx, max_val = self.get_peak()
            tmp_idx = self.current_idx - buffer_size + max_idx
            self.window_locked = append(self.window_locked, self.current_idx)
            if tmp_idx == self.last_mhat:
                self.state = st.locked
            elif abs(tmp_idx - (self.last_mhat + frame_size)) <= 1:
                self.state = st.locked
                self.llast_mhat = self.last_mhat
                self.last_mhat = tmp_idx
                self.mhats = append(self.mhats, self.last_mhat)
            elif abs(tmp_idx - (self.llast_mhat + frame_size)) <= 1:
                self.state = st.locked
                self.last_mhat = tmp_idx
                self.mhats = append(self.mhats, self.last_mhat)
            else:
                self.state = st.oneMissed
                self.weird = tmp_idx
        elif self.state == st.oneMissed:
            if self.current_idx < self.weird + frame_size:
                self.state == st.oneMissed
            else:
                avg_val, max_idx, max_val = self.get_peak()
                tmp_idx = self.current_idx - buffer_size + max_idx
                if abs(tmp_idx - (self.last_mhat + frame_size)) <= 1:
                    self.state = st.locked
                    self.llast_mhat = self.last_mhat
                    self.last_weird = None
                    self.last_mhat = tmp_idx
                    self.mhats = append(self.mhats, self.last_mhat)
                    self.weird = None
                else:
                    self.state = st.twoMissed
                    self.weird = tmp_idx
        elif self.state == st.twoMissed:
            if self.current_idx < self.weird + frame_size:
                self.state == st.twoMissed
            else:
                avg_val, max_idx, max_val = self.get_peak()
                tmp_idx = self.current_idx - buffer_size + max_idx
                if (
                    abs(tmp_idx - (self.last_mhat + 2 * frame_size)) <= 1
                    or abs(tmp_idx - (self.last_mhat + frame_size)) <= 1
                ):
                    self.state = st.locked
                    self.llast_mhat = self.last_mhat
                    self.last_weird = None
                    self.last_mhat = tmp_idx
                    self.mhats = append(self.mhats, self.last_mhat)
                    self.weird = None
                else:
                    self.state = st.onePeak
                    self.last_weird = self.last_mhat
                    self.last_mhat = tmp_idx
                    self.mhats = append(self.mhats, self.last_mhat)
                    self.llast_mhat = None
                    self.weird = None
                    self.lock_end = self.current_idx
                    lock_pair = array([self.lock_begin,self.lock_end])
                    self.lock_pairs = append(self.lock_pairs, lock_pair)
        else:
            print("Error state")

def plot_window(sm, samples, peak, num_plots=3):

    # plot the surrounding pointd of the peak (256 samples on either side)
    data_to_correlate = samples[peak-256:peak+513]
    possible_preamble = zeros(513)
    x = range(512)
    for index in x:
        possible_preamble[index] = L0_single(data_to_correlate[index:index+256])
    x = x + peak - 257
    max_y = max(possible_preamble)
    max_x = argmax(possible_preamble) + peak - 257

    # plot the window of correlated data the state machine is seeing
    window = sm.correlations
    win_x = range(sm.current_idx - buffer_size, sm.current_idx - 255)
    winMax_y = max(window)
    winMax_x = argmax((window)) + sm.current_idx - buffer_size

    # actually do the plotting
    fig = figure(figsize=(12,9), dpi=80)
    axs = fig.subplots(num_plots,1)
    axs[0].plot(x, possible_preamble[:512])
    axs[0].text(max_x-64, max_y + 500, "max = (" + str(max_x) + "," + str(max_y) + ")")
    axs[0].grid(True)
    if num_plots > 1:
        axs[1].plot(win_x, window)
        axs[1].text(winMax_x-128, max_y + 500, "max = (" + str(winMax_x) + "," + str(winMax_y) + ")")
        axs[1].axvspan(winMax_x-256, winMax_x+256, color='red', alpha=0.5)
        axs[1].grid(True)
        if num_plots > 2:
            complete_idx = range(-256,sm.current_idx-255)
            axs[2].plot(complete_idx, sm.all_correlated_samples)
            axs[2].grid(True)
            # if prod(sm.lock_pairs.shape) > 0:
            #     for i in range(shape(sm.lock_pairs)[0]):
            #         axs[2].axvspan(sm.lock_pairs[i][0], sm.lock_pairs[i][1], color='yellow', alpha=.2)
            for num in complete_idx:
                if num % 6656 == 4504:
                    is_mhat = False
                    for mhat in sm.mhats:
                        if (abs(mhat - num) <= 2):
                            is_mhat = True
                    if is_mhat:
                        axs[2].axvspan(num-256, num+256, color='green', alpha=0.4)
                    else:
                        axs[2].axvspan(num-256, num+256, color='red', alpha=0.4)
                    
    show()  

def main():
    sm = FSSM()
    samples = genfromtxt("../data/9_Aug_Wired.csv", dtype=complex).flatten()
    # samples = genfromtxt("../MATLAB/data_a_aug25.csv", dtype=complex).flatten() 
    last_state = st.init
    count = 0
    last_mhat = None
    weird = None
    for sample in samples:
        # Print every 1000 samples so we know where it's at
        if not count % 1000:
            print(count)
        count += 1

        # Print the state and wait
        if sm.state != last_state:
            last_state = sm.state
            print(sm.state, end="\n")

            # input()

        # Tick
        sm.tick(sample)

        # Print mhats and weirds
        if sm.last_mhat != last_mhat:
            last_mhat = sm.last_mhat
            print(f"last_mhat: {last_mhat}")
            plot_window(sm, samples, last_mhat)

        if sm.weird != weird:
            weird = sm.weird
            if weird is not None:
                print(f"weird: {weird}")
                plot_window(sm, samples, weird)

                



if __name__ == "__main__":
    main()