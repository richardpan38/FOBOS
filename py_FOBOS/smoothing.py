import numpy as np
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import pylab as pl
import sys
import scipy.ndimage.filters as ndimage

def border_plot():
        #Generating the border for the inner Mini-IFU
        border_row = []
        border_column = []

        for i in range(1 , 10):
                for j in range(1 , 10):
                        if 75 <= (i ** 2) + (j ** 2) <= (radi ** 2):
                                border_row.append(i)
                                border_column.append(j)
        for i in range(len(circle_row)):
            border_row.append( border_row[i])
            border_column.append(-border_column[i])
        for j in range(len(border_row)):
                border_row.append(-border_row[j])
                border_column.append(border_column[j])

        #Generating the borders for the outer Mini-IFUs

        outer_border_row = []
        outer_border_column = []
        for i in range(1 , 26):
                for j in range(1 , 26):
                        if 600 < (i ** 2) + (j ** 2) <= (radi_outer ** 2):
                                outer_border_row.append(i)
                                outer_border_column.append(j)
        for i in range(9 , 23):
                for j in range(1 , 23):
                        if -.5 < i - (j* np.sqrt(3)) < .5:
                                outer_border_row.append(i)
                                outer_border_column.append(j)


        for i in range(len(outer_border_row)):
            outer_border_row.append( outer_circle_row[i])
            outer_border_column.append(-outer_circle_column[i])
        for j in range(len(outer_border_row)):
                outer_border_row.append(-outer_border_row[j])
                outer_border_column.append(outer_border_column[j])

        for i in range(0,26):
                if (radi ** 2) < (i ** 2) <= (radi_outer ** 2):
                        outer_border_row.append(0)
                        outer_border_column.append(i)
                        outer_border_row.append(0)
                        outer_border_column.append(-i)


        #Dividing the outer border to respective fibers
        border_fiber1_x = []
        border_fiber1_y = []
        border_fiber2_x = []
        border_fiber2_y = []
        border_fiber3_x = []
        border_fiber3_y = []
        border_fiber4_x = []
        border_fiber4_y = []
        border_fiber5_x = []
        border_fiber5_y = []
        border_fiber6_x = []
        border_fiber6_y = []

        for i in range(0 , len(outer_border_row)):
                if outer_border_row[i] > 0 and outer_border_column[i] > 0:
                        if 0  <  (outer_border_row[i]/outer_border_column[i])  < 1.65:
                                border_fiber1_x.append(outer_border_column[i])
                                border_fiber1_y.append(outer_border_row[i])
                        else:
                                border_fiber2_x.append(outer_border_column[i])
                                border_fiber2_y.append(outer_border_row[i])
                if outer_border_row[i] < 0 and outer_border_column[i] > 0:
                        if 0  <  (-outer_border_row[i]/outer_border_column[i])  <= 1.85:
                                border_fiber6_x.append(outer_border_column[i])
                                border_fiber6_y.append(outer_border_row[i])
                        else:
                                border_fiber5_x.append(outer_border_column[i])
                                border_fiber5_y.append(outer_border_row[i])
                if outer_border_row[i] > 0 and outer_border_column[i] < 0:
                        if 0  <  (outer_border_row[i]/-outer_border_column[i])  <= 1.85:
                                border_fiber3_x.append(outer_border_column[i])
                                border_fiber3_y.append(outer_border_row[i])
                        else:
                                border_fiber2_x.append(outer_border_column[i])
                                border_fiber2_y.append(outer_border_row[i])
                if outer_border_row[i] < 0 and outer_border_column[i] < 0:
                        if 0  <  (-outer_border_row[i]/-outer_border_column[i])  < 1.65:
                                border_fiber4_x.append(outer_border_column[i])
                                border_fiber4_y.append(outer_border_row[i])
                        else:
                                fiber5_x.append(outer_border_column[i])
                                fiber5_y.append(outer_border_row[i])
        for i in range(0 , len(outer_border_row)):
                if outer_border_row[i] == 0 and outer_border_column[i] > 0:
                        border_fiber1_x.append(outer_border_column[i])
                        border_fiber1_y.append(outer_border_row[i])
                if outer_border_row[i] > 0 and outer_border_column[i] == 0:
                        border_fiber2_x.append(outer_border_column[i])
                        border_fiber2_y.append(outer_border_row[i])
                if outer_border_row[i] == 0 and outer_border_column[i] < 0:
                        border_fiber4_x.append(outer_border_column[i])
                        border_fiber4_y.append(outer_border_row[i])
                if outer_border_row[i] < 0 and outer_border_column[i] == 0:
                        border_fiber5_x.append(outer_border_column[i])
                        border_fiber5_y.append(outer_border_row[i])

        #Creating fibers for each object
        global border_fiber0_xcoords
        global border_fiber0_ycoords
        global border_fiber1_xcoords
        global border_fiber1_ycoords
        global border_fiber2_xcoords
        global border_fiber2_ycoords
        global border_fiber3_xcoords
        global border_fiber3_ycoords
        global border_fiber4_xcoords
        global border_fiber4_ycoords
        global border_fiber5_xcoords
        global border_fiber5_ycoords
        global border_fiber6_xcoords
        global border_fiber6_ycoords
        border_fiber1_xcoords = []
        border_fiber1_ycoords = []
        border_fiber2_xcoords = []
        border_fiber2_ycoords = []
        border_fiber3_xcoords = []
        border_fiber3_ycoords = []
        border_fiber4_xcoords = []
        border_fiber4_ycoords = []
        border_fiber5_xcoords = []
        border_fiber5_ycoords = []
        border_fiber6_xcoords = []
        border_fiber6_ycoords = []
        border_fiber1_row = [50 + y for y in border_fiber1_y]
        border_fiber1_column = [50 + x for x in border_fiber1_x]
        border_fiber2_row = [50 + y for y in border_fiber2_y]
        border_fiber2_column = [50 + x for x in border_fiber2_x]
        border_fiber3_row = [50 + y for y in border_fiber3_y]
        border_fiber3_column = [50 + x for x in border_fiber3_x]
        border_fiber4_row = [50 + y for y in border_fiber4_y]
        border_fiber4_column = [50 + x for x in border_fiber4_x]
        border_fiber5_row = [50 + y for y in border_fiber5_y]
        border_fiber5_column = [50 + x for x in border_fiber5_x]
        border_fiber6_row = [50 + y for y in border_fiber6_y]
        border_fiber6_column = [50 + x for x in border_fiber6_x]  
        border_fiber1_xcoords.extend(border_fiber1_column)
        border_fiber1_ycoords.extend(border_fiber1_row)
        border_fiber2_xcoords.extend(border_fiber2_column)
        border_fiber2_ycoords.extend(border_fiber2_row)
        border_fiber3_xcoords.extend(border_fiber3_column)
        border_fiber3_ycoords.extend(border_fiber3_row)
        border_fiber4_xcoords.extend(border_fiber4_column)
        border_fiber4_ycoords.extend(border_fiber4_row)
        border_fiber5_xcoords.extend(border_fiber5_column)
        border_fiber5_ycoords.extend(border_fiber5_row)
        border_fiber6_xcoords.extend(border_fiber6_column)
        border_fiber6_ycoords.extend(border_fiber6_row)

        #Generating all the x,y coordinates for each object

        border_fiber0_xcoords = []
        border_fiber0_ycoords = []
        fiber0_border_row = [ y + 50  for y in border_row]
        fiber0_border_column = [ x + 50 for x in border_column]
        border_fiber0_xcoords.append(x_coords[target_index])
        border_fiber0_ycoords.append(y_coords[target_index])
        border_fiber0_xcoords.extend(fiber0_border_column)
        border_fiber0_ycoords.extend(fiber0_border_row)

def smoothing(x):
       FWHM = x
       sigma1 = fwhm2sigma(FWHM)/A_P
       flux0_gaussian = scidata[(y_coords[target_index] - 50):(y_coords[target_index] + 50), (x_coords[target_index] - 50): (y_coords[target_index] + 50)] 
       gaussian_0 = ndimage.gaussian_filter(flux0_gaussian, sigma = sigma1)
       plt.scatter(border_fiber0_xcoords, border_fiber0_ycoords)
       plt.scatter(border_fiber1_xcoords, border_fiber1_ycoords)
       plt.scatter(border_fiber2_xcoords, border_fiber2_ycoords)
       plt.scatter(border_fiber3_xcoords, border_fiber3_ycoords)
       plt.scatter(border_fiber4_xcoords, border_fiber4_ycoords)
       plt.scatter(border_fiber5_xcoords, border_fiber5_ycoords)
       plt.scatter(border_fiber6_xcoords, border_fiber6_ycoords)
       plt.imshow(gaussian_0, cmap = 'gray')
       plt.colorbar()
       plt.show()



