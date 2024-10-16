import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def ramped_reg_wave(T, d_t, H, T_wave, T_pos_sprung, T_width_sprung):
    """creates ramped regular wave as array and Figure, smothing with 3-wise convolution
    input:
        T: float, windwo size
        d_t: float, timestep
        H: foat, waveheiht
        T_wave: float, Waveperiad
        T_pos_sprung: float, Position middel of rampig
        T_width_sprung: float, Width of ramp

    output:
        elev: array, elevation of ramped wave
        fig: Figure, Figure of elevation of ramped wave and convolution
    """

    N_sprung = int(T_pos_sprung / d_t)
    N_width_sprung = int(T_width_sprung / d_t)
    t = np.arange(0, T, d_t)
    N = len(t)

    sprung = np.concatenate(((np.zeros(N_sprung)), np.ones(N - N_sprung)))

    box_flatten = np.ones(N_width_sprung) / N_width_sprung

    sprung = np.convolve(sprung, box_flatten, mode='same')
    sprung = np.convolve(sprung, box_flatten, mode='same')
    sprung_flatten = np.convolve(sprung, box_flatten, mode='same')
    sprung_flatten[N_sprung + int(N_width_sprung):-1] = 1

    elev = H / 2 * np.sin(1 / T_wave * 2 * np.pi * t)
    elev_flatten = elev * sprung_flatten

    fig, ax = plt.subplots(1, 1)

    ax.plot(t, sprung_flatten)
    ax.plot(t, elev_flatten)
    ax.axvline(T_pos_sprung, color="blue")
    ax.axvline(T_pos_sprung - T_width_sprung / 2, color="blue", linestyle="dashed")
    ax.axvline(T_pos_sprung + T_width_sprung / 2, color="blue", linestyle="dashed")
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Waveelevation [m]')
    ax.set_title(f'Ramped Wave with: H={H}, T_wave={T_wave}, T_width={T_width_sprung}, Tpos_sprung={T_pos_sprung}')

    t = t.reshape(N, 1)
    elev_flatten = elev_flatten.reshape(N, 1)

    elev = np.concatenate((t, elev_flatten), axis=1)

    return elev, fig


def calculate_damp(x, threshold=None, bothSides=False, decay=True):
    """ Detect maxima in a signal, computes the log deg based on it

    input:
        x: Seires, Signal of Occilation

    optional:
        threshold: float or None, Threshold for detecting maxima, std max/3
        bothSides: bool, If False, only maxima are detected
        decay: bool, If True, logaritmic decay is asumed

    output:
        t_ref: array: referenced time of logdec, T, sigma, D
        median_logdec: float, Median logaritmic decrement
        median_T: float, Median Waveperiod of Occilation
        median_sigma: float, Median decay constant of Occilation
        median_D: float, Median of damping ratio of Occilation
        logdec: array, logaritmic decrement between peaks
        T: array, Waveperiod between peaks
        sigma: array, decay constant between peaks
        D: array, damping ratio between peaks

    Theory:
    Logarithmic decrement:
        delta = 1/N log [ x(t) / x(t+N T_d)]  = 2 pi zeta / sqrt(1-zeta^2)

    Damping ratio:
            zeta = delta / sqrt( 4 pi^2 + delta^2 )

    Damped period, frequency:
        Td = 2pi / omega_d
        omegad = omega_0 sqrt(1-zeta**2)

    Damping exponent:
        alpha = zeta omega_0 = delta/ T_d
    """

    def indexes(y, thres=0.3, min_dist=1, thres_abs=False):
        """Peak detection routine.

        Finds the numeric index of the peaks in *y* by taking its first order difference. By using
        *thres* and *min_dist* parameters, it is possible to reduce the number of
        detected peaks. *y* must be signed.

        Parameters
        ----------
        y : ndarray (signed)
            1D amplitude data to search for peaks.
        thres : float, defining threshold. Only the peaks with amplitude higher than the
            threshold will be detected.
            if thres_abs is False: between [0., 1.], normalized threshold.
        min_dist : int
            Minimum distance between each detected peak. The peak with the highest
            amplitude is preferred to satisfy this constraint.
        thres_abs: boolean
            If True, the thres value will be interpreted as an absolute value, instead of
            a normalized threshold.

        Returns
        -------
        ndarray
            Array containing the numeric indexes of the peaks that were detected
        """
        if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
            raise ValueError("y must be signed")

        if not thres_abs:
            thres = thres * (np.max(y) - np.min(y)) + np.min(y)

        min_dist = int(min_dist)

        # compute first order difference
        dy = np.diff(y)

        # propagate left and right values successively to fill all plateau pixels (0-value)
        zeros, = np.where(dy == 0)

        # check if the signal is totally flat
        if len(zeros) == len(y) - 1:
            return np.array([])

        if len(zeros):
            # compute first order difference of zero indexes
            zeros_diff = np.diff(zeros)
            # check when zeros are not chained together
            zeros_diff_not_one, = np.add(np.where(zeros_diff != 1), 1)
            # make an array of the chained zero indexes
            zero_plateaus = np.split(zeros, zeros_diff_not_one)

            # fix if leftmost value in dy is zero
            if zero_plateaus[0][0] == 0:
                dy[zero_plateaus[0]] = dy[zero_plateaus[0][-1] + 1]
                zero_plateaus.pop(0)

            # fix if rightmost value of dy is zero
            if len(zero_plateaus) and zero_plateaus[-1][-1] == len(dy) - 1:
                dy[zero_plateaus[-1]] = dy[zero_plateaus[-1][0] - 1]
                zero_plateaus.pop(-1)

            # for each chain of zero indexes
            for plateau in zero_plateaus:
                median = np.median(plateau)
                # set leftmost values to leftmost non zero values
                dy[plateau[plateau < median]] = dy[plateau[0] - 1]
                # set rightmost and middle values to rightmost non zero values
                dy[plateau[plateau >= median]] = dy[plateau[-1] + 1]

        # find the peaks by using the first order difference
        peaks = np.where((np.hstack([dy, 0.]) < 0.)
                         & (np.hstack([0., dy]) > 0.)
                         & (np.greater(y, thres)))[0]

        # handle multiple peaks, respecting the minimum distance
        if peaks.size > 1 and min_dist > 1:
            highest = peaks[np.argsort(y[peaks])][::-1]
            rem = np.ones(y.size, dtype=bool)
            rem[peaks] = False

            for peak in highest:
                if not rem[peak]:
                    sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                    rem[sl] = True
                    rem[peak] = False

            peaks = np.arange(y.size)[~rem]

        return peaks

    values = x.values
    values_mean = np.mean(values)
    values = values - values_mean

    t = x.index

    if bothSides:
        ldPos, iTPos, stdPos, IPos, vldPos = calculate_damp(x, threshold=threshold, decay=decay)
        ldNeg, iTNeg, stdNeg, INeg, vldNeg = calculate_damp(-x, threshold=threshold, decay=decay)
        return (ldPos + ldNeg) / 2, (iTPos + iTNeg) / 2, (stdPos + stdNeg) / 2, (IPos, INeg), (vldPos, vldNeg)

    if threshold is None:
        threshold = np.mean(abs(x - np.mean(x))) / 3

    I = indexes(values, thres=threshold, min_dist=1, thres_abs=True)
    t_ref = (t[I[1:]]+t[I[:-1]])/2

    T = t[I[1:]]-t[I[:-1]]

    median_T = np.median(T)

    vn = np.arange(0, len(I) - 1) + 1

    # Quick And Dirty Way using one as ref and assuming all periods were found
    if decay:
        # For a decay we take the first peak as a reference
        logdec = 1 / vn * np.log(values[I[0]] / values[I[1:]])  # Logarithmic decrement
    else:
        # For negative damping we take the last peak as a reference
        logdec = 1 / vn * np.log(values[I[-2::-1]] / values[I[-1]])  # Logarithmic decrement

    median_logdec = np.mean(logdec)
    sigma = logdec/T
    median_sigma = np.median(sigma)

    D = logdec/(np.pi*2)
    median_D= median_logdec/(np.pi*2)

    return t_ref, median_logdec, median_T, median_sigma, median_D, logdec, T, sigma, D

