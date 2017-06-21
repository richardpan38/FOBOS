import numpy as np
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import pylab as pl
import sys
import scipy.ndimage.filters as ndimage

hdulist = fits.open('/Users/RichardP/research/FOBOS/Flux/Samples/acs_I_030mas_088_sci.fits')


#In degrees!
RA = hdulist[0].header['RA_Targ']
DEC = hdulist[0].header['DEC_Targ']

#Arcsecond:Pixel Conversion factor and Flux Area
#.57 arcseconds is .015833333 degrees
Dia = .57
radi_arcs = Dia/2
outer_radi =.01583333333 * np.sqrt(7)
A_P = .03
radi = ((Dia/A_P)/2)
fiber_area  = np.pi * (((Dia/A_P)/2)**2)
xy_max = radi/(np.sqrt(2))

#Each pixel
scidata = hdulist[0].data

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

catalog = fits.open('../../Flux/Samples/photoz_vers2.0_010312_UltraVISTA2016.fits')
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


objects_index = np.where(np.logical_and(catalog_RA  <= RA_max , catalog_RA >= RA_min))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_DEC <= DEC_max , catalog_DEC >= DEC_min)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Imag > 0 ,catalog_Umag > 0)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Zmag > 0 ,catalog_Vmag > 0)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Rmag > 0 ,catalog_Gmag > 0)))
objects_index = np.intersect1d(objects_index, np.where(catalog_Bmag > 0))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_zgal > 0 ,catalog_zilb > 0)))
objects_index = np.intersect1d(objects_index, np.where(np.logical_and(catalog_Bmag > 0 ,catalog_chis > 0)))
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
objects2 = box_4.wcs_pix2world(objects,1)
x_coords = [row[0] for row in objects]
y_coords = [row[1] for row in objects]


#Creating the first inner fiber
# square_row = []
# square_column = []

# for i in range(1 , 50):
#         for j in range(1 , 50):
#             square_row.append(i)
#             square_column.append(j)
# for i in range(len(square_row)):
#     square_row.append(square_row[i])
#     square_column.append(-square_column[i])
# for j in range(len(square_row)):
#         square_row.append(-square_row[j])
#         square_column.append(square_column[j])
# for i in range(1,50):
#     square_row.append(i)
#     square_column.append(0)
#     square_row.append(0)
#     square_column.append(i)
#     square_row.append(-i)
#     square_column.append(0)
#     square_row.append(0)
#     square_column.append(-i)
# square_row.append(0)
# square_column.append(0)

circle_row = []
circle_column = []

for i in range(1 , 10):
        for j in range(1 , 10):
                if 75 <= (i ** 2) + (j ** 2) <= (radi ** 2):
                        circle_row.append(i)
                        circle_column.append(j)
for i in range(len(circle_row)):
    circle_row.append( circle_row[i])
    circle_column.append(-circle_column[i])
for j in range(len(circle_row)):
        circle_row.append(-circle_row[j])
        circle_column.append(circle_column[j])



# for i in range(1,10):
#         if 0 <= (i ** 2) <= (radi ** 2):
#                 circle_row.append(i)
#                 circle_column.append(0)
#                 circle_row.append(0)
#                 circle_column.append(i)
#                 circle_row.append(-i)
#                 circle_column.append(0)
#                 circle_row.append(0)
#                 circle_column.append(-i)
# circle_row.append(0)
# circle_column.append(0)

radi_outer = radi * np.sqrt(7)
outer_circle_row = []
outer_circle_column = []
for i in range(1 , 26):
        for j in range(1 , 26):
                if 600 < (i ** 2) + (j ** 2) <= (radi_outer ** 2):
                        outer_circle_row.append(i)
                        outer_circle_column.append(j)
for i in range(9 , 23):
        for j in range(1 , 23):
                if -.5 < i - (j* np.sqrt(3)) < .5:
                        outer_circle_row.append(i)
                        outer_circle_column.append(j)
                

for i in range(len(outer_circle_row)):
    outer_circle_row.append( outer_circle_row[i])
    outer_circle_column.append(-outer_circle_column[i])
for j in range(len(outer_circle_row)):
        outer_circle_row.append(-outer_circle_row[j])
        outer_circle_column.append(outer_circle_column[j])

for i in range(0,26):
        if (radi ** 2) < (i ** 2) <= (radi_outer ** 2):
#                outer_circle_row.append(i)
#               outer_circle_column.append(0)
                outer_circle_row.append(0)
                outer_circle_column.append(i)
#               outer_circle_row.append(-i)
#               outer_circle_column.append(0)
                outer_circle_row.append(0)
                outer_circle_column.append(-i)
fiber1_x = []
fiber1_y = []
fiber2_x = []
fiber2_y = []
fiber3_x = []
fiber3_y = []
fiber4_x = []
fiber4_y = []
fiber5_x = []
fiber5_y = []
fiber6_x = []
fiber6_y = []

for i in range(0 , len(outer_circle_row)):
        if outer_circle_row[i] > 0 and outer_circle_column[i] > 0:
                if 0  <  (outer_circle_row[i]/outer_circle_column[i])  < 1.65:
                        fiber1_x.append(outer_circle_column[i])
                        fiber1_y.append(outer_circle_row[i])
                else:
                        fiber2_x.append(outer_circle_column[i])
                        fiber2_y.append(outer_circle_row[i])
        if outer_circle_row[i] < 0 and outer_circle_column[i] > 0:
                if 0  <  (-outer_circle_row[i]/outer_circle_column[i])  <= 1.85:
                        fiber6_x.append(outer_circle_column[i])
                        fiber6_y.append(outer_circle_row[i])
                else:
                        fiber5_x.append(outer_circle_column[i])
                        fiber5_y.append(outer_circle_row[i])
        if outer_circle_row[i] > 0 and outer_circle_column[i] < 0:
                if 0  <  (outer_circle_row[i]/-outer_circle_column[i])  <= 1.85:
                        fiber3_x.append(outer_circle_column[i])
                        fiber3_y.append(outer_circle_row[i])
                else:
                        fiber2_x.append(outer_circle_column[i])
                        fiber2_y.append(outer_circle_row[i])
        if outer_circle_row[i] < 0 and outer_circle_column[i] < 0:
                if 0  <  (-outer_circle_row[i]/-outer_circle_column[i])  < 1.65:
                        fiber4_x.append(outer_circle_column[i])
                        fiber4_y.append(outer_circle_row[i])
                else:
                        fiber5_x.append(outer_circle_column[i])
                        fiber5_y.append(outer_circle_row[i])
for i in range(0 , len(outer_circle_row)):
        if outer_circle_row[i] == 0 and outer_circle_column[i] > 0:
                fiber1_x.append(outer_circle_column[i])
                fiber1_y.append(outer_circle_row[i])
        if outer_circle_row[i] > 0 and outer_circle_column[i] == 0:
                fiber2_x.append(outer_circle_column[i])
                fiber2_y.append(outer_circle_row[i])
        if outer_circle_row[i] == 0 and outer_circle_column[i] < 0:
                fiber4_x.append(outer_circle_column[i])
                fiber4_y.append(outer_circle_row[i])
        if outer_circle_row[i] < 0 and outer_circle_column[i] == 0:
                fiber5_x.append(outer_circle_column[i])
                fiber5_y.append(outer_circle_row[i])

#Creating fibers for each object

fiber1_xcoords = []
fiber1_ycoords = []
fiber2_xcoords = []
fiber2_ycoords = []
fiber3_xcoords = []
fiber3_ycoords = []
fiber4_xcoords = []
fiber4_ycoords = []
fiber5_xcoords = []
fiber5_ycoords = []
fiber6_xcoords = []
fiber6_ycoords = []

object_number = 9

# outer_fiber1_row = [y_coords[object_number] + y for y in fiber1_y]
# outer_fiber1_column = [x_coords[object_number] + x for x in fiber1_x]
# outer_fiber2_row = [y_coords[object_number] + y for y in fiber2_y]
# outer_fiber2_column = [x_coords[object_number] + x for x in fiber2_x]
# outer_fiber3_row = [y_coords[object_number] + y for y in fiber3_y]
# outer_fiber3_column = [x_coords[object_number] + x for x in fiber3_x]
# outer_fiber4_row = [y_coords[object_number] + y for y in fiber4_y]
# outer_fiber4_column = [x_coords[object_number] + x for x in fiber4_x]
# outer_fiber5_row = [y_coords[object_number] + y for y in fiber5_y]
# outer_fiber5_column = [x_coords[object_number] + x for x in fiber5_x]
# outer_fiber6_row = [y_coords[object_number] + y for y in fiber6_y]
# outer_fiber6_column = [x_coords[object_number] + x for x in fiber6_x]
outer_fiber1_row = [50 + y for y in fiber1_y]
outer_fiber1_column = [50 + x for x in fiber1_x]
outer_fiber2_row = [50 + y for y in fiber2_y]
outer_fiber2_column = [50 + x for x in fiber2_x]
outer_fiber3_row = [50 + y for y in fiber3_y]
outer_fiber3_column = [50 + x for x in fiber3_x]
outer_fiber4_row = [50 + y for y in fiber4_y]
outer_fiber4_column = [50 + x for x in fiber4_x]
outer_fiber5_row = [50 + y for y in fiber5_y]
outer_fiber5_column = [50 + x for x in fiber5_x]
outer_fiber6_row = [50 + y for y in fiber6_y]
outer_fiber6_column = [50 + x for x in fiber6_x]  
fiber1_xcoords.extend(outer_fiber1_column)
fiber1_ycoords.extend(outer_fiber1_row)
fiber2_xcoords.extend(outer_fiber2_column)
fiber2_ycoords.extend(outer_fiber2_row)
fiber3_xcoords.extend(outer_fiber3_column)
fiber3_ycoords.extend(outer_fiber3_row)
fiber4_xcoords.extend(outer_fiber4_column)
fiber4_ycoords.extend(outer_fiber4_row)
fiber5_xcoords.extend(outer_fiber5_column)
fiber5_ycoords.extend(outer_fiber5_row)
fiber6_xcoords.extend(outer_fiber6_column)
fiber6_ycoords.extend(outer_fiber6_row)
# plt.scatter(fiber1_x, fiber1_y)
# plt.scatter(fiber2_x, fiber2_y)
# plt.scatter(fiber3_x, fiber3_y)
# plt.scatter(fiber4_x, fiber4_y)
# plt.scatter(fiber5_x, fiber5_y)
# plt.scatter(fiber6_x, fiber6_y)
# plt.scatter(circle_column, circle_row)
# plt.show()
# sys.exit()



#Generating all the x,y coordinates for each object

object_x_coords = []
object_y_coords = []
fiber0_xcoords = []
fiber0_ycoords = []

fiber_row = [ y + 50  for y in circle_row]
fiber_column = [ x + 50 for x in circle_column]
fiber0_xcoords.append(x_coords[object_number])
fiber0_ycoords.append(y_coords[object_number])
fiber0_xcoords.extend(fiber_column)
fiber0_ycoords.extend(fiber_row)
plt.scatter(fiber0_xcoords, fiber0_ycoords)
plt.scatter(fiber1_xcoords, fiber1_ycoords)
plt.scatter(fiber2_xcoords, fiber2_ycoords)
plt.scatter(fiber3_xcoords, fiber3_ycoords)
plt.scatter(fiber4_xcoords, fiber4_ycoords)
plt.scatter(fiber5_xcoords, fiber5_ycoords)
plt.scatter(fiber6_xcoords, fiber6_ycoords)
#for object_number in range(0, len(x_coords)):
#        fiber_row = [y_coords[object_number] + y for y in square_row]
#        fiber_column = [x_coords[object_number] + x for x in square_column]
#        fiber0_xcoords.append(x_coords[object_number])
#        fiber0_ycoords.append(y_coords[object_number])
#        fiber0_xcoords.extend(fiber_column)
#        fiber0_ycoords.extend(fiber_row)
        
#8733
#15388
flux0_gaussian = hdulist[0].data[15338:15438 , 8683:8783]


#Gaussian Smoothing for first point
def fwhm2sigma(fwhm):
        return fwhm / np.sqrt(8 * np.log(2))
FWHM = .5
sigma1 = fwhm2sigma(FWHM)/A_P
gaussian_0 = ndimage.gaussian_filter(flux0_gaussian, sigma = sigma1)
plt.imshow(gaussian_0, cmap = 'gray')
plt.colorbar()
plt.show()

FWHM = .75
sigma1 = fwhm2sigma(FWHM)/A_P
gaussian_0 = ndimage.gaussian_filter(flux0_gaussian, sigma = sigma1)
plt.scatter(fiber0_xcoords, fiber0_ycoords)
plt.scatter(fiber1_xcoords, fiber1_ycoords)
plt.scatter(fiber2_xcoords, fiber2_ycoords)
plt.scatter(fiber3_xcoords, fiber3_ycoords)
plt.scatter(fiber4_xcoords, fiber4_ycoords)
plt.scatter(fiber5_xcoords, fiber5_ycoords)
plt.scatter(fiber6_xcoords, fiber6_ycoords)
plt.imshow(gaussian_0, cmap = 'gray')
plt.colorbar()
plt.show()

FWHM = 1
sigma1 = fwhm2sigma(FWHM)/A_P
gaussian_0 = ndimage.gaussian_filter(flux0_gaussian, sigma = sigma1)
plt.scatter(fiber0_xcoords, fiber0_ycoords)
plt.scatter(fiber1_xcoords, fiber1_ycoords)
plt.scatter(fiber2_xcoords, fiber2_ycoords)
plt.scatter(fiber3_xcoords, fiber3_ycoords)
plt.scatter(fiber4_xcoords, fiber4_ycoords)
plt.scatter(fiber5_xcoords, fiber5_ycoords)
plt.scatter(fiber6_xcoords, fiber6_ycoords)
plt.imshow(gaussian_0, cmap = 'gray')
plt.colorbar()
plt.show()

FWHM = 1.25
sigma1 = fwhm2sigma(FWHM)/A_P
gaussian_0 = ndimage.gaussian_filter(flux0_gaussian, sigma = sigma1)
plt.scatter(fiber0_xcoords, fiber0_ycoords)
plt.scatter(fiber1_xcoords, fiber1_ycoords)
plt.scatter(fiber2_xcoords, fiber2_ycoords)
plt.scatter(fiber3_xcoords, fiber3_ycoords)
plt.scatter(fiber4_xcoords, fiber4_ycoords)
plt.scatter(fiber5_xcoords, fiber5_ycoords)
plt.scatter(fiber6_xcoords, fiber6_ycoords)
plt.imshow(gaussian_0, cmap = 'gray')
plt.colorbar()
plt.show()

FWHM = 1.5
sigma1 = fwhm2sigma(FWHM)/A_P
gaussian_0 = ndimage.gaussian_filter(flux0_gaussian, sigma = sigma1)
plt.scatter(fiber0_xcoords, fiber0_ycoords)
plt.scatter(fiber1_xcoords, fiber1_ycoords)
plt.scatter(fiber2_xcoords, fiber2_ycoords)
plt.scatter(fiber3_xcoords, fiber3_ycoords)
plt.scatter(fiber4_xcoords, fiber4_ycoords)
plt.scatter(fiber5_xcoords, fiber5_ycoords)
plt.scatter(fiber6_xcoords, fiber6_ycoords)
plt.imshow(gaussian_0, cmap = 'gray')
plt.colorbar()
plt.show()

#Conversions to ABMAG
def ABMAG_Convert(FluxVals):
    ABMAG_list = []
    for i in range(0, len(FluxVals)):
        temp = abs(FluxVals[i]) /4
        ABMAG =  ((-2.5 * np.log(temp)) +23.9)
        ABMAG_list.append(ABMAG)
        


#Plotting just a couple objects
test_x = []

#test_x = hdulist[0].data[15290:15544 , 6700:6954].tolist()
# print(np.shape(test_x))

#print(x_coords[0])
#print(y_coords[0])
#plt.plot(test_gaussian)
#plt.colorbar()
#plt.scatter(plot_y, objects_Imag)
# plt.plot(objects_Imag)
#plt.show()



