"""Adapted from Matlab code: https://terpconnect.umd.edu/~toh/spectrum/findpeaksL.m

Original Matlab version by:
    T. C. O'Haver, 1995, 2014

Translated to Python in 2017 by Paul Freeman
"""
from itertools import combinations
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class Peak:
    """Represents a peak"""
    def __init__(self, x_data, y_data, center, amplitude, num_points):
        center = max(x_data[num_points // 2], min(center, x_data[len(x_data) - num_points // 2 - 2]))

        center_index = (np.abs(x_data - center)).argmin()
        self.start = max(center_index - num_points // 2, 0)
        self.end = min(self.start + num_points, len(x_data) - 1)

        peak_x_values = x_data[self.start:self.end]
        peak_y_values = y_data[self.start:self.end]
        if len(peak_x_values) != num_points:
            raise ValueError('peak_x_values is short: center was {}'.format(center))
        if len(peak_y_values) != num_points:
            raise ValueError('peak_y_values is short: center was {}'.format(center))

        # use poly fit to make guesses
        a, b, c = tuple(np.polyfit(peak_x_values, np.reciprocal(peak_y_values), 2)[:3])
        center_guess = -b / (2*a)
        half_width_guess = np.sqrt(abs((4*a*c-b**2) / a) / np.sqrt(abs(a))) / 2
        if x_data[num_points // 2] > center_guess or center_guess > x_data[len(x_data) - num_points // 2 - 2]:
            center_guess = x_data[center_index]
            half_width_guess = 1

        # fit curve to Lorentzian function
        center_min = peak_x_values[0]
        center_max = peak_x_values[-1]
        width_min = 3 * np.abs(x_data[-1] - x_data[0]) / len(x_data)
        width_max = 0.125 * np.abs(x_data[-1] - x_data[0])
        amp_min = amplitude * 0.75
        amp_max = amplitude / 0.75
        guess = max(center_min, min(center_guess, center_max)), max(width_min, min(half_width_guess, width_max)), max(amp_min, min(amplitude, amp_max))
        bounds = ([center_min, width_min, amp_min],
                  [center_max, width_max, amp_max])
        try:
            popt, pcov = curve_fit(
                _loren,
                peak_x_values,
                peak_y_values,
                p0=guess,
                bounds=bounds)
        except ValueError:
            print('guess = {}'.format(guess))
            print('bounds = {}'.format(bounds))
            raise
        self.center, half_width, self.amplitude = tuple(popt)
        self.width = half_width * 2
        self.err = np.sqrt(np.diag(pcov))

    def __lt__(self, other):
        return self.center < other.center

class PeakPicker:
    """Interactive plot to pick peaks"""
    def __init__(self, x_data, y_data, peak_width):
        self.x_data = x_data
        self.y_data = y_data
        self.peak_width = max(3, peak_width)
        self.peaks = []
        self.visible_peaks = True
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.scatter(self.x_data, self.y_data, s=3)
        self.ax.set_title('Peak picker')
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
        #plt.show()

    def onkey(self, event):
        if event.key == 'tab':
            self.visible_peaks = not self.visible_peaks
            self.draw_peaks()

    def onclick(self, event):
        if event.inaxes != self.ax:
            return
        if event.button == 1:
            try:
                self.peaks.append(Peak(self.x_data, self.y_data, event.xdata, event.ydata, self.peak_width))
            except RuntimeError as err:
                print('Could not find peak: {}'.format(err))
                return
        elif event.button == 3:
            x_values = np.array([peak.center for peak in self.peaks])
            dists = np.abs(x_values - event.xdata)
            self.peaks.pop(dists.argmin())
        self.draw_peaks()

    def toolsbar(self):
	toolbar = plt.get_current_fig_manager().toolbar
	if toolbar.mode!='':
		return True
	return False

    def draw_peaks(self):
        print('{} peaks:'.format(len(self.peaks)))
        self.peaks.sort()
        self.ax.clear()
        if self.visible_peaks:
            self.ax.scatter(self.x_data, self.y_data, s=3)
        self.refine2(5)
        #self.refine3(2)
        fit = np.zeros(len(self.x_data))
        for peak in self.peaks:
            print('Center: {:.3f} (error {:.3f}), Width: {:.3f} (error {:.3f}), Amplitude: {:.3f} (error {:.3f})'.format(
                peak.center, peak.err[0], peak.width, peak.err[1], peak.amplitude, peak.err[2]))
            curve = loren(self.x_data, peak.center, peak.width, peak.amplitude)
            fit += curve
            if self.visible_peaks:
                self.ax.plot(self.x_data, curve, 'orange')
                self.ax.plot((peak.center - peak.width / 2, peak.center + peak.width / 2),
                             (peak.amplitude / 2, peak.amplitude / 2),
                             'green')
                self.ax.annotate(s='width {:.1f}'.format(peak.width),
                                 xy=(peak.center, peak.amplitude / 2),
                                 xytext=(10, -20),
                                 textcoords="offset points",
                                 bbox=dict(boxstyle='square, pad=0.5', fc='white', alpha=0.7)
                                 )
                self.ax.plot((peak.center, peak.center),
                             (0, peak.amplitude),
                             'green', picker=5)
                self.ax.annotate(s='peak ({:.1f}, {:.1f})'.format(peak.center, peak.amplitude),
                                 xy=(peak.center, peak.amplitude),
                                 xytext=(10, 20),
                                 textcoords="offset points",
                                 bbox=dict(boxstyle='square, pad=0.5', fc='white', alpha=0.7)
                                 )
        ssd = np.sum((fit - self.y_data)**2)
        print('Total error (SSD): {:.6f}'.format(ssd))
        if self.visible_peaks:
            self.ax.plot(self.x_data, fit, 'red')
        else:
            self.ax.scatter(self.x_data, self.y_data - fit, s=3)
        self.fig.canvas.draw()

    def refine2(self, passes):
        for _ in range(passes):
            for peak1, peak2 in combinations(self.peaks, 2):
                #start = None
                #end = None
                #if peak1 is not self.peaks[0]:
                start = peak1.start
                #if peak2 is not self.peaks[-1]:
                end = peak2.end
                peaks_x_values = self.x_data[start:end]
                isolate_peaks = self.y_data[start:end].copy()
                for peak in self.peaks:
                    if peak is not peak1 and peak is not peak2:
                        isolate_peaks -= loren(peaks_x_values, peak.center, peak.width, peak.amplitude)
                peaks_y_values = isolate_peaks

                # refit curves
                width_min = 3 * np.abs(self.x_data[-1] - self.x_data[0]) / len(self.x_data)
                width_max = 0.125 * np.abs(self.x_data[-1] - self.x_data[0])
                lower_bound1 = [self.x_data[peak1.start], width_min, 0.9 * peak1.amplitude]
                lower_bound2 = [self.x_data[peak2.start], width_min, 0.9 * peak2.amplitude]
                upper_bound1 = [self.x_data[peak1.end], width_max, 1.1 * peak1.amplitude]
                upper_bound2 = [self.x_data[peak2.end], width_max, 1.1 * peak2.amplitude]
                guess = (peak1.center, peak1.width / 2, peak1.amplitude,
                         peak2.center, peak2.width / 2, peak2.amplitude)
                try:
                    popt, pcov = curve_fit(
                        lambda x, ctr1, wid1, amp1, ctr2, wid2, amp2: (
                            _loren(x, ctr1, wid1, amp1)
                            + _loren(x, ctr2, wid2, amp2)),
                        peaks_x_values,
                        peaks_y_values,
                        p0=guess,
                        bounds=(lower_bound1 + lower_bound2, upper_bound1 + upper_bound2)
                        )
                except ValueError:
                    return
                except RuntimeError:
                    return
                (peak1.center, half_width1, peak1.amplitude,
                 peak2.center, half_width2, peak2.amplitude) = tuple(popt)
                center_index1 = (np.abs(self.x_data - peak1.center)).argmin()
                #peak1.start = max(center_index1 - self.peak_width // 2, 0)
                #peak1.end = min(peak1.start + self.peak_width, len(self.x_data) - 1)
                center_index2 = (np.abs(self.x_data - peak2.center)).argmin()
                #peak2.start = max(center_index2 - self.peak_width // 2, 0)
                #peak2.end = min(peak2.start + self.peak_width, len(self.x_data) - 1)
                peak1.width = half_width1 * 2
                peak2.width = half_width2 * 2
                peak1.err = np.sqrt(np.diag(pcov))[:3]
                peak2.err = np.sqrt(np.diag(pcov))[3:]

    def refine3(self, passes):
        for _ in range(passes):
            for peak1, peak2, peak3 in combinations(self.peaks, 3):
                #start = None
                #end = None
                #if peak1 is not self.peaks[0]:
                start = peak1.start
                #if peak3 is not self.peaks[-1]:
                end = peak2.end
                peaks_x_values = self.x_data[start:end]
                peaks_y_values = self.y_data[start:end].copy()
                for peak in self.peaks:
                    if peak is not peak1 and peak is not peak2 and peak is not peak3:
                        peaks_y_values -= loren(peaks_x_values, peak.center, peak.width, peak.amplitude)

                # refit curves
                width_min = 3 * np.abs(self.x_data[-1] - self.x_data[0]) / len(self.x_data)
                width_max = 0.125 * np.abs(self.x_data[-1] - self.x_data[0])
                lower_bound1 = [self.x_data[peak1.start], width_min, 0.9 * peak1.amplitude]
                lower_bound2 = [self.x_data[peak2.start], width_min, 0.9 * peak2.amplitude]
                lower_bound3 = [self.x_data[peak3.start], width_min, 0.9 * peak3.amplitude]
                upper_bound1 = [self.x_data[peak1.end], width_max, 1.1 * peak1.amplitude]
                upper_bound2 = [self.x_data[peak2.end], width_max, 1.1 * peak2.amplitude]
                upper_bound3 = [self.x_data[peak3.end], width_max, 1.1 * peak3.amplitude]
                guess = (peak1.center, peak1.width / 2, peak1.amplitude,
                         peak2.center, peak2.width / 2, peak2.amplitude,
                         peak3.center, peak3.width / 2, peak3.amplitude)
                try:
                    popt, pcov = curve_fit(
                        lambda x, ctr1, wid1, amp1, ctr2, wid2, amp2, ctr3, wid3, amp3: (
                            _loren(x, ctr1, wid1, amp1)
                            + _loren(x, ctr2, wid2, amp2)
                            + _loren(x, ctr3, wid3, amp3)),
                        peaks_x_values,
                        peaks_y_values,
                        p0=guess,
                        bounds=(lower_bound1 + lower_bound2 + lower_bound3,
                                upper_bound1 + upper_bound2 + upper_bound3)
                        )
                except ValueError:
                    return
                except RuntimeError:
                    return
                (peak1.center, half_width1, peak1.amplitude,
                 peak2.center, half_width2, peak2.amplitude,
                 peak3.center, half_width3, peak3.amplitude) = tuple(popt)
                center_index1 = (np.abs(self.x_data - peak1.center)).argmin()
                #peak1.start = max(center_index1 - self.peak_width // 2, 0)
                #peak1.end = min(peak1.start + self.peak_width, len(self.x_data) - 1)
                center_index2 = (np.abs(self.x_data - peak2.center)).argmin()
                #peak2.start = max(center_index2 - self.peak_width // 2, 0)
                #peak2.end = min(peak2.start + self.peak_width, len(self.x_data) - 1)
                center_index3 = (np.abs(self.x_data - peak3.center)).argmin()
                #peak3.start = max(center_index3 - self.peak_width // 2, 0)
                #peak3.end = min(peak3.start + self.peak_width, len(self.x_data) - 1)
                peak1.width = half_width1 * 2
                peak2.width = half_width2 * 2
                peak3.width = half_width3 * 2
                peak1.err = np.sqrt(np.diag(pcov))[:3]
                peak2.err = np.sqrt(np.diag(pcov))[3:6]
                peak3.err = np.sqrt(np.diag(pcov))[6:]

def loren(x_values, center, width, amp=1):
    """Calculate a Lorentzian curve"""
    half_width = width / 2.0
    return np.array([_loren(x, center, half_width, amp) for x in x_values])

def _loren(x, center, half_width, amp=1):
    """Calculate a Lorentzian value"""
    return amp / (1 + ((x - center) / half_width)**2)


if __name__ == '__main__':
    X_DATA = np.arange(1, 100, 0.002)
    Y_DATA = (loren(X_DATA, 20, 4.1, 1.3) +
              loren(X_DATA, 25, 5.0, 1.0) +
              loren(X_DATA, 40, 5.8, 1.2) +
              loren(X_DATA, 50, 7.1, 2.1) +
              loren(X_DATA, 80, 3.9, 0.9) +
              loren(X_DATA, 83, 8.7, 0.4))
    NOISE = 0.3 * (np.random.rand(len(X_DATA)) - 0.5)
    print('Test')
    print('Curve #1: Center: 20, Width: 4.1, Amplitude: 1.3')
    print('Curve #2: Center: 25, Width: 5.0, Amplitude: 1.0')
    print('Curve #3: Center: 40, Width: 5.8, Amplitude: 1.2')
    print('Curve #4: Center: 50, Width: 7.1, Amplitude: 2.1')
    print('Curve #5: Center: 80, Width: 3.9, Amplitude: 0.9')
    print('Curve #6: Center: 83, Width: 1.7, Amplitude: 0.6')
    print('Noise: +/-0.1')
    p = PeakPicker(X_DATA, Y_DATA + NOISE, 1700)
    #plt.show()
