from __future__ import division, print_function

import numpy
from astropy import wcs
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy.constants as sp

filename = 'Users/RichardP/research/FOBOS/Flux/Samples/acs_I_030mas_088_sci.fits'

hdulist = fits.open('/Users/RichardP/research/FOBOS/Flux/Samples/acs_I_030mas_088_sci.fits')
lambda1 = hdulist[0].header["photplam"]
#image_data = hdulist[0].data
#plt.imshow(image_data, cmap='gray')
#plt.show()


outer_radi =.01583333333 * np.sqrt(7)
#Each pixel
scidata = hdulist[0].data
#photflam = hdulist[0].header['photflam']
#exptime = hdulist[0].header['exptime']
#scidata *= photflam / exptime


#Instead convert the 4 points of the box to create the dimensions

box_4= wcs.WCS(hdulist[0].header)
box_coords = np.array([(25,25),(25,20454),(20455, 25),(20455,20455)])
Box_RA_DEC = box_4.wcs_pix2world(box_coords, 1)
RA_m = Box_RA_DEC[2,0].tolist()
RA_ma = Box_RA_DEC[0,0].tolist()
DEC_m = Box_RA_DEC[0,1].tolist()
DEC_ma = Box_RA_DEC[1,1].tolist()
RA_min = RA_m + outer_radi
RA_max = RA_ma - outer_radi
DEC_min = DEC_m + outer_radi
DEC_max = DEC_ma - outer_radi

#From the catalog find all the targets within the range of the image
#Generate a list of of all the RA and DEC
#Using if statement look through each ID to find values in the range
#of the RA and DEC of the specific .fits file
#Finding only the relevant catalog_RA and catalog_DEC

catalog = fits.open('../Samples/photoz_vers2.0_010312_UltraVISTA2016.fits')
catalog_wcs = wcs.WCS(catalog[1].header)
catalog_id = catalog[1].data['id']
catalog_RA = np.array(catalog[1].data['RA'])
catalog_DEC = np.array(catalog[1].data['DEC'])
catalog_Imag = np.array(catalog[1].data['imag'])
catalog_Umag = np.array(catalog[1].data['umag'])
catalog_Zmag = np.array(catalog[1].data['zmag'])
catalog_Rmag = np.array(catalog[1].data['rmag'])
catalog_Vmag = np.array(catalog[1].data['vmag'])
catalog_Gmag = np.array(catalog[1].data['gmag'])
catalog_Bmag = np.array(catalog[1].data['bmag'])
catalog_zgal = np.array(catalog[1].data['zgal'])
catalog_zilb = np.array(catalog[1].data['z_ilb'])
catalog_chis = np.array(catalog[1].data['chi_sec'])
catalog_zsec = np.array(catalog[1].data['zsec'])

# objects_index = []

# for i in range(0, len(catalog_RA)):
#         if RA_min <= catalog_RA[i] <= RA_max and DEC_min <= catalog_DEC[i] <= DEC_max and catalog_Imag[i] > 0 and catalog_Umag[i] > 0 and catalog_Zmag[i] > 0 and catalog_Vmag[i] > 0 and catalog_Rmag[i] > 0 and catalog_Gmag[i] > 0 and catalog_Bmag[i] > 0:
#                 objects_index.append(i)
# print(len(objects_index))
# sys.exit()

objects_index = np.where(np.logical_and(catalog_RA  <= RA_max , catalog_RA >= RA_min))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_DEC <= DEC_max , catalog_DEC >= DEC_min)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Imag > 0 ,catalog_Umag > 0)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Zmag > 0 ,catalog_Vmag > 0)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Rmag > 0 ,catalog_Gmag > 0)))
objects_index = np.intersect1d(objects_index, np.where(catalog_Bmag > 0))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_zgal > 0 ,catalog_zilb > 0)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Bmag > 0 ,catalog_chis > 0)))
#objects_index = np.intersect1d(objects_index, np.where(catalog_zsec > 0))
objects_id = catalog_id[objects_index]
objects_RA = catalog_RA[objects_index]
objects_DEC = catalog_DEC[objects_index]
objects_Imag = catalog_Imag[objects_index]
objects_Umag = catalog_Umag[objects_index]
objects_Zmag = catalog_Zmag[objects_index]
objects_Vmag = catalog_Vmag[objects_index]
objects_Rmag = catalog_Rmag[objects_index]
objects_Gmag = catalog_Gmag[objects_index]
objects_Bmag =  catalog_Bmag[objects_index]


#Convert catalog WCS to PIXELS MUCH MORE EFFICIENT
#Most x,y coords were in decimals so I made them into integers. Some of them rounded upwards
x = np.array([(objects_RA[0], objects_DEC[0])], np.float)
for i in range(1, len(objects_index)):
        y = np.array([(objects_RA[i], objects_DEC[i])])
        x = np.vstack(( x , y))
objects = (box_4.wcs_world2pix(x , 1))
objects = objects.astype(int)
#print(objects[0])
objects2 = box_4.wcs_pix2world(objects,1)
#print(objects2[0])
x_coords = [row[0] for row in objects]
y_coords = [row[1] for row in objects]
#print(catalog[1].data[273726])

#Plot Flux versus Imag to see relationship
# DS9 Plotting Tools



     
#Finding the points that fall within the inner circle fiber
#This is for any arbitrary circle of radius (radi).
circle_row = []
circle_column = []


#Creating the first inner fiber

A_P = .03
radi = 1.5/.03

circle_row = []
circle_column = []

for i in range(1 , 50):
        for j in range(1 , 50):
                if 0 <= (i ** 2) + (j ** 2) <= (radi ** 2):
                        circle_row.append(i)
                        circle_column.append(j)
for i in range(len(circle_row)):
    circle_row.append( circle_row[i])
    circle_column.append(-circle_column[i])
for j in range(len(circle_row)):
        circle_row.append(-circle_row[j])
        circle_column.append(circle_column[j])



for i in range(1,50):
        if 0 <= (i ** 2) <= (radi ** 2):
                circle_row.append(i)
                circle_column.append(0)
                circle_row.append(0)
                circle_column.append(i)
                circle_row.append(-i)
                circle_column.append(0)
                circle_row.append(0)
                circle_column.append(-i)
circle_row.append(0)
circle_column.append(0)

#Generating all the x,y coordinates for each object

object_x_coords = []
object_y_coords = []
fiber0_xcoords = []
fiber0_ycoords = []

for object_number in range(0, len(x_coords)):
        fiber_row = [y_coords[object_number] + y for y in circle_row]
        fiber_column = [x_coords[object_number] + x for x in circle_column]
        fiber0_xcoords.append(x_coords[object_number])
        fiber0_ycoords.append(y_coords[object_number])
        fiber0_xcoords.extend(fiber_column)
        fiber0_ycoords.extend(fiber_row)
        
#NOTE: EACH OBJECT has 293 pixels for fiber0
#Generate flux values for each object's fiber0

flux0_values = []
for i in range(0,len(fiber0_xcoords)):
        flux0_values.append(scidata[fiber0_ycoords[i] , fiber0_xcoords[i]])

fiber_sum = {}

for i in range(0, len(x_coords)):
        counter0 = i * 7841

        fiber_sum[str(objects_id[i]) + '_fiber0'.format(i)] = sum(flux0_values[counter0 + j] for j in range(7841))


#Created an array formatted as so to sum flux of each object. 
object_info = np.array([(0, 0)])

for i in range(0, len(x_coords)):
        flux_sum = float(fiber_sum[str(objects_id[i]) + '_fiber0'])
        index = np.array([(objects_id[i], flux_sum)])
        object_info = np.vstack((object_info , index))

plot_y = [float(row[1]) for row in object_info]
del plot_y[0]
print(len(plot_y))
#plt.scatter(plot_y, objects_Imag)
#plt.show()
#ABMAG = -2.5 Log Fν - 48.60. The Fv = flux * photoflam / exptime.
photflam = hdulist[0].header['photflam']
exptime = hdulist[0].header['exptime']
print(photflam)
print(exptime)
#print(len(plot_y))
#plot = plot_y[1] * 5
#print(len(plot))
#plot_y *= photflam/exptime
objects_Imag = objects_Imag.tolist()
plot_imag = []
for i in range(0, len(plot_y)):
        temp = abs(plot_y[i]) /4
        ABMAG =  ((-2.5 * np.log(temp)) +23.9)
#        Flux = (10 ** ((23.9 - objects_Imag[i])/2.5)) * 4.75
        plot_imag.append(ABMAG)
#        plot_imag.append(Flux)
#print((plot_imag))
plt.scatter(objects_Imag,  plot_imag)
#axes = plt.gca()
#axes.set_xlim([15,30])
#axes.set_ylim([15,30])
plt.show()





