import numpy as np
from matplotlib import pyplot as plt

data = np.load("../../../temp_data/tests_data.npy")
geo_locs = np.load("../../../temp_data/tests_geo_loc.npy")
lbls = np.load("../../../temp_data/tests_labels.npy")
dt = 0.0005

tt = np.arange(0, data.shape[1]) * dt
for wig, loc in zip(data, geo_locs):
    if loc[1] == 0 and int(2 * loc[0]) % 2 == 0 and (int(2 * loc[0]) // 2) % 5 == 0:
        plt.plot(wig * 2000 + loc[0], tt, label=loc[0], c="k")

plt.gca().invert_yaxis()
plt.xlabel("x offset (meters)")
plt.ylabel("time (sec)")
plt.title("simulated data example - 10 sources, on y = 0")
plt.show()
