import numpy as np
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import pylab as pl
import sys
import scipy.ndimage.filters as ndimage

def ABMAG_Convert(FluxVals):
    ABMAG_list = []
    for i in range(0, len(FluxVals)):
        temp = abs(FluxVals[i]) /4
        ABMAG =  ((-2.5 * np.log(temp)) +23.9)
        ABMAG_list.append(ABMAG)
        print(ABMag_list)
