import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics.traveltime import TravelTimeManager
from matplotlib.widgets import Cursor


import pygimli.meshtools as mt
from pygimli.physics import traveltime
from pygimli.viewer.pv import drawSensors
from pygimli.physics.seismics import solvePressureWave

from pygimli.physics import seismics as seis


pyvista = pg.optImport("pyvista")

pg.utils.units.quants['vel']['cMap'] = 'inferno_r'

SEED=42
np.random.seed(SEED)


def create_mesh():
    pass


def main():
    """

    :return:
    """
    create_mesh()


if __name__ == "__main__":
    main()
