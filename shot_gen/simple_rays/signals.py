import numpy as np


def ricker(t, t0=0, s=1, a=1):
    """
    create a ricker wavelet at a certain time, centered around t0, with a chosen spread
    :param t: the time at which we're calculating the wavelet (can be an numpy array of time)
    :param t0: the time the wavelet occurs at
    :param s: the spread of the wavelet
    :param a: the amplitude of the wavelet
    :return: the corresponding ricker wavelet
    """
    t = (t - t0) / s
    amp = 3 * np.sqrt(np.pi) * s
    amp = a * 2 / np.sqrt(amp)

    f_exp = - (t ** 2) / 2

    res = t ** 2
    res = amp * (1 - res)
    return res * np.exp(f_exp)


def ricker_tr(s=1):
    """
    calculate the radius of the time window the wavelet takes place in
    :param s: the spread of the wavelet
    :return: the radius of the time window in which the wavelet appears
    """
    return 5.5 * s


def sinc(t, t0=0, s=1, a=1):
    """ricker
    create a sinc signal at a certain time, centered around t0, with a chosen spread
    :param t: the time at which we're calculating the wavelet (can be an numpy array of time)
    :param t0: the time the signal occurs at
    :param s: the spread of the wavelet
    :param a: the amplitude of the wavelet
    :return: the corresponding sinc signal
    """
    t = (t - t0) / s
    return a * np.sinc(t)


def sinc_tr(s=1):
    """
    calculate the radius of the time window the wavelet takes place in
    :param s: the spread of the wavelet
    :return: the radius of the time window in which the wavelet appears
    """
    return 50 * s


def vlad(t, t0=0, s=1, a=1):
    """
    create a vlad wavelet at a certain time, centered around t0, with a chosen spread
    :param t: the time at which we're calculating the wavelet (can be an numpy array of time)
    :param t0: the time the wavelet occurs at
    :param s: the spread of the wavelet
    :param a: the amplitude of the wavelet
    :return: the corresponding vlad wavelet
    """
    t = (t - t0) / s
    amp = np.abs(1 * t)
    amp = np.minimum(amp, np.ones_like(amp))
    amp = a * ((1 - amp) ** 2)
    f = 2 * np.pi * t
    return amp * np.sin(f)


def vlad_tr(s=1):
    """
    calculate the radius of the time window the wavelet takes place in
    :param s: the spread of the wavelet
    :return: the radius of the time window in which the wavelet appears
    """
    return 1.0 * s
