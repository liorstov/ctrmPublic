import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import argparse
import os

SEED=42
np.random.seed(SEED)


def gen_a_shot(f, t0, params, geophs, loc, t_span, vp=None, vs=None):
    """
    given a geophoes geometry, a disturbance location and type,
    will generate a a timeline with only this disturbance
    :param f: a function describing the disturbance, of form f(t, *params)
    :param t0: the time the disturbance occurs at
    :param params: a list of params to pass to f() after the time
    :param geophs: a numpy array containing the locations of the geophones
    :param loc: the location of the disturbance
    :param t_span: the amount of time the recorded in the shot
    :param vp: the speed of the pressure wave. If None, will be ignored.
    :param vs: the speed of the sheer wave. If None, will be ignored.
    :return: a numpy array containing the shot data of the disturbance
    """
    # make sure that we're working with some kind of wave
    assert (vp is not None) or (vs is not None), "can't have both vp and vs being None"
    # save the number of geophones in memory
    n_geophs = geophs.shape[0]
    # find the distance of each geophone from the source
    dists = np.sqrt(((np.array([[loc[0], loc[1], loc[2]]] * n_geophs) - geophs)**2).sum(axis=1))
    # calculate the decay of the amplitude by the time it arrives to each geophone
    cos_phi = np.abs(geophs[:, 2] - loc[2]) / dists
    amps = cos_phi / (dists**2)
    # TODO: find max window size, after which the disturbance has decayed (simplifies computation)
    # TODO: add sheer waves support
    shots = np.zeros((dists.shape[0], t_span[1] - t_span[0]))
    for v in [vp, vs]:
        if v is not None:
            shots += gen_shot_for_speed(f, v, dists, t_span, t0, amps, params)
    return shots


def gen_shot_for_speed(f, v, dists, t_span, t0, amps, params):
    """
    generate shot data for 1 speed mode
    :param f: a function describing the disturbance, of form f(t, *params)
    :param v: the speed of the waves in the medium
    :param dists: the distances of the wave source from each sensor
    :param t_span: the time span the experiment takes place in: [t_0, t_f]
    :param t0: the time the disturbance occurs at
    :param amps: the amplitude decay of the the signal to the sensor
    :param params: additional parameters the source function takes
    :return: shot data for each sensor as a numpy matrix
    """
    # find the time the disturbance takes to travel to each geophone
    dt = dists / v
    # create an object to hold the output data
    shots = np.zeros((dists.shape[0], t_span[1] - t_span[0]))
    for t in np.arange(t_span[0], t_span[1]):
        shots[:, t] = amps * f(t - dt, t0, *params)
    return shots


def gen_shots(f, param_range, geophs, loc_range, N, max_subs, break_range, sub_range, vp=None, vs=None, seed=42):
    """
    generate a sequence of disturbances and record shot data
    :param f: a function describing the disturbance, of form f(t, t0, *params)
    :param param_range: 2 lists describing the range of params that can be passed to to f()
    after the time variable: [params_0, params_f]
    :param geophs: a numpy array containing the locations of the geophones
    :param loc_range: 2 lists describing the bounds of the possible disturbances:
    [(x_0, y_0, z_0), (x_f, y_f, z_f)]
    :param N: the number of disturbances
    :param max_subs: how many disturbances at most can make each disturbance
    :param break_range: a range of time to wait between disturbances
    :param sub_range: a range of time to wait between sub disturbances
    :param vp: speed of pressure waves (if None, will be ignored)
    :param vs: speed of sheer waves (if None, will be ignored)
    :param seed: used to set the randomness of the model
    :return: a numpy array containing the shot data of the disturbance
    """
    # calculate the moment each disturbance takes place, and total sim time
    # TODO: implement seed usage, or delete it as an input
    T = 0
    dis_times = []
    dis_dt = []
    n_geophs = geophs.shape[0]
    for _ in range(N):
        T += np.random.randint(*break_range)
        n_subs = np.random.randint(1, max_subs + 1)
        t_wait = np.random.randint(*sub_range)
        dis_dt.append(t_wait)
        sub_dists = []
        for _ in range(n_subs):
            T += t_wait
            sub_dists.append(T)
        dis_times.append(sub_dists)
    T += np.random.randint(*break_range)

    # create object to hold output data
    outp = np.zeros((n_geophs, T))
    y_res = np.zeros((4, T))

    # create a random disturbance each dis moment
    for disturb_range, dis_wait in tqdm(zip(dis_times, dis_dt)):
        params = []
        for a, b in zip(param_range[0], param_range[1]):
            params.append(a + ((b - a) * np.random.random()))
        loc = []
        for a, b in zip(loc_range[0], loc_range[1]):
            loc.append(a + ((b - a) * np.random.random()))

        tt0 = disturb_range[0]
        ttf = disturb_range[-1] + dis_wait
        y_res[:, tt0:ttf] = np.array([[1] + loc] * (ttf - tt0)).T

        for disturb_t in disturb_range:
            outp += gen_a_shot(f, disturb_t, params, geophs, loc, [0, T], vp=vp, vs=vs)
    return outp, y_res


def create_geophone_grid(x0, xf, dx, y0, yf, dy, z=0):
    """
    create a grid of geophones
    :param x0: the start of the x axis
    :param xf: the end of the x axis
    :param dx: the step size of the x axis
    :param y0: the start of the y axis
    :param yf: the end of the y axis
    :param dy: the step size of the y axis
    :param z: the height of the grid (either a num, or a list)
    :return: the grid as a numpy array (as a numpy array)
    """
    geophs = []
    try:
        _ = iter(z)
    except TypeError:
        z = [z]
    for x in np.arange(x0, xf + (0.5 * dx), dx):
        for y in np.arange(y0, yf + (0.5 * dy), dy):
            for iz in z:
                geophs.append(np.array([x, y, iz]))
    geophs = np.array(geophs)
    return geophs


def create_geophone_square(s0, sf, ds, z=0):
    """
    create a square grid of geophones (x0=y0, xf=yf)
    :param s0: the start of the grid
    :param sf: the end of the grid
    :param ds: step size on the grid
    :param z: the height of the grid
    :return: the grid as a numpy array (as a numpy array)
    """
    return create_geophone_grid(s0, sf, ds, s0, sf, ds, z)


def get_args():
    """
    :return: argument parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--s0", help="start of the grid", type=float)
    parser.add_argument("--sf", help="end of the grid", type=float)
    parser.add_argument("--ds", help="grid resolution", type=float)
    parser.add_argument("-vp", help="speed of pressure waves", type=float, default=None)
    parser.add_argument("-vs", help="speed of sheer waves", type=float, default=None)
    parser.add_argument("-N", help="number of disturbances", type=int)
    parser.add_argument("-z", help="a list of z values for the geophone grid", type=float, default=0.0, nargs='+')
    parser.add_argument("--z0", help="start of disturbance depth range", type=float, default=0.0)
    parser.add_argument("--zf", help="end of disturbance depth range", type=float, default=-50.0)
    parser.add_argument("--max-subs", help="max disturbances per disturbance", type=int, default=1)

    parser.add_argument("--break-range0", help="time to wait between disturbances - start", type=int, default=400)
    parser.add_argument("--break-rangef", help="time to wait between disturbances - end", type=int, default=500)
    parser.add_argument("--sub-range0", help="time to wait between sub disturbances - start", type=int, default=50)
    parser.add_argument("--sub-rangef", help="time to wait between  subdisturbances - end", type=int, default=100)
    parser.add_argument("--SEED", help="a seed for randomness", type=int, default=SEED)

    parser.add_argument("--A0", help="range of amplitudes - start", type=float, default=0.9)
    parser.add_argument("--Af", help="range of amplitudes - end", type=float, default=1.1)
    parser.add_argument("--w0", help="range of wave size - start", type=float, default=0.1)
    parser.add_argument("--wf", help="range of wave size - end", type=float, default=0.3)

    parser.add_argument('--show-data', help="show graph of example data", action='store_true')

    return parser


def main():
    """
    generate shot data, along with a label file for location and time
    :return: 1 is successful
    """
    args = get_args().parse_args()
    if SEED != args.SEED:
        np.random.seed(args.SEED)

    data_file_name = os.path.abspath("./data.npy")
    lbls_file_name = os.path.abspath("./labels.npy")
    geos_file_name = os.path.abspath("./geo_loc.npy")
    s0 = args.s0
    sf = args.sf
    ds = args.ds
    geophs = create_geophone_square(s0, sf, ds, z=args.z)
    sig_sinc = lambda t, t0, A, w: A * np.sinc(w * (t - t0))
    data, y_true = gen_shots(sig_sinc,
                             [[args.A0, args.w0], [args.Af, args.wf]],
                             geophs,
                             [(s0, s0, args.z0), (sf, sf, args.zf)],
                             N=args.N,
                             max_subs=args.max_subs,
                             break_range=[args.break_range0, args.break_rangef],
                             sub_range=[args.sub_range0, args.sub_rangef],
                             vp=args.vp,
                             vs=args.vs,
                             seed=args.SEED)
    if args.show_data:
        print(f"the shape of data is {data.shape}")
        print(f"the shape of the results is {y_true.shape}")
        plt.plot(np.arange(0, data.shape[1]), data[-1, :] / (data[-1, :].max()), label="last shot")
        plt.plot(np.arange(0, data.shape[1]), data[0, :] / (data[0, :].max()), label="first shot")
        plt.plot(np.arange(0, data.shape[1]), y_true[0, :], label="times of signal")
        plt.xlabel("time")
        plt.ylabel("normalized shots")
        plt.title("example shots - normalized")
        plt.legend()
        plt.show()
    np.save(data_file_name, data)
    np.save(lbls_file_name, y_true)
    np.save(geos_file_name, geophs)
    return 1


if __name__ == "__main__":
    main()
