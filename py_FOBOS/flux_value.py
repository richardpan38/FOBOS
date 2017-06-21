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

catalog = fits.open('/Users/RichardP/Research/FOBOS/Samples/photoz_vers2.0_010312_UltraVISTA2016.fits')
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
objects2 = box_4.wcs_pix2world(objects,1)
x_coords = [row[0] for row in objects]
y_coords = [row[1] for row in objects]

#Creating the first inner fiber
circle_row = []
circle_column = []

for i in range(1 , 10):
        for j in range(1 , 10):
                if 0 <= (i ** 2) + (j ** 2) <= (radi ** 2):
                        circle_row.append(i)
                        circle_column.append(j)
for i in range(len(circle_row)):
    circle_row.append( circle_row[i])
    circle_column.append(-circle_column[i])
for j in range(len(circle_row)):
        circle_row.append(-circle_row[j])
        circle_column.append(circle_column[j])

for i in range(1,10):
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
#Radius of the outer circle
#All the points on the outer circle

radi_outer = radi * np.sqrt(7)
outer_circle_row = []
outer_circle_column = []
for i in range(1 , 26):
        for j in range(1 , 26):
                if (radi ** 2)  < (i ** 2) + (j ** 2) <= (radi_outer ** 2):
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
                outer_circle_row.append(i)
                outer_circle_column.append(0)
                outer_circle_row.append(0)
                outer_circle_column.append(i)
                outer_circle_row.append(-i)
                outer_circle_column.append(0)
                outer_circle_row.append(0)
                outer_circle_column.append(-i)


#Splitting up each fiber before adding object coordinates
#Remade the bounds to produce bounds closer to 283:283
                
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


for object_number in range(0 , len(y_coords)):
        outer_fiber1_row = [y_coords[object_number] + y for y in fiber1_y]
        outer_fiber1_column = [x_coords[object_number] + x for x in fiber1_x]
        outer_fiber2_row = [y_coords[object_number] + y for y in fiber2_y]
        outer_fiber2_column = [x_coords[object_number] + x for x in fiber2_x]
        outer_fiber3_row = [y_coords[object_number] + y for y in fiber3_y]
        outer_fiber3_column = [x_coords[object_number] + x for x in fiber3_x]
        outer_fiber4_row = [y_coords[object_number] + y for y in fiber4_y]
        outer_fiber4_column = [x_coords[object_number] + x for x in fiber4_x]
        outer_fiber5_row = [y_coords[object_number] + y for y in fiber5_y]
        outer_fiber5_column = [x_coords[object_number] + x for x in fiber5_x]
        outer_fiber6_row = [y_coords[object_number] + y for y in fiber6_y]
        outer_fiber6_column = [x_coords[object_number] + x for x in fiber6_x]   
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


#Given pixel coordinates I set new variables for each individual pixel's flux value
#Then I will create a new list of all the objects with variables set to their fiber's total flux
#Fiber1/4: 285
#Fiber2/5: 282
#Fiber3/6: 283

flux1_values = []
flux2_values = []
flux3_values = []
flux4_values = []
flux5_values = []
flux6_values = []

for i in range(0,len(fiber1_xcoords)):
        flux1_values.append(scidata[fiber1_ycoords[i] , fiber1_xcoords[i]])
for i in range(0,len(fiber2_xcoords)):
        flux2_values.append(scidata[fiber2_ycoords[i] , fiber2_xcoords[i]])
for i in range(0,len(fiber3_xcoords)):
        flux3_values.append(scidata[fiber3_ycoords[i] , fiber3_xcoords[i]])
for i in range(0,len(fiber4_xcoords)):
        flux4_values.append(scidata[fiber4_ycoords[i] , fiber4_xcoords[i]])
for i in range(0,len(fiber5_xcoords)):
        flux5_values.append(scidata[fiber5_ycoords[i] , fiber5_xcoords[i]])
for i in range(0,len(fiber6_xcoords)):
        flux6_values.append(scidata[fiber6_ycoords[i] , fiber6_xcoords[i]])

fiber_sum = {}
        
for i in range(0, len(x_coords)):
        counter0 = i * 293
        counter1 = i * 285
        counter2 = i * 282
        counter3 = i * 283
        fiber_sum[str(objects_id[i]) + '_fiber0'.format(i)] = sum(flux0_values[counter0 + j] for j in range(293))
        fiber_sum[str(objects_id[i]) + '_fiber1'.format(i)] = sum(flux1_values[counter1 + k] for k in range(285))
        fiber_sum[str(objects_id[i]) + '_fiber2'.format(i)] = sum(flux2_values[counter2 + l] for l in range(282))
        fiber_sum[str(objects_id[i]) + '_fiber3'.format(i)] = sum(flux3_values[counter3 + m] for m in range(283))
        fiber_sum[str(objects_id[i]) + '_fiber4'.format(i)] = sum(flux4_values[counter1 + n] for n in range(285))
        fiber_sum[str(objects_id[i]) + '_fiber5'.format(i)] = sum(flux5_values[counter2 + o] for o in range(282))
        fiber_sum[str(objects_id[i]) + '_fiber6'.format(i)] = sum(flux6_values[counter3 + p] for p in range(283))


#Created an array formatted as so to sum flux of each object. 
object_info = np.array([(0 , 0)])

for i in range(0, len(x_coords)):
        flux_sum = [float(fiber_sum[str(objects_id[i]) + '_fiber0']) , float(fiber_sum[str(objects_id[i]) + '_fiber1']), float(fiber_sum[str(objects_id[i]) + '_fiber2']) , float(fiber_sum[str(objects_id[i]) + '_fiber3']) , float(fiber_sum[str(objects_id[i]) + '_fiber4']) , float(fiber_sum[str(objects_id[i]) + '_fiber5']) , float(fiber_sum[str(objects_id[i]) + '_fiber6'])]
        index = np.array([(objects_id[i], sum(flux_sum))])
        object_info = np.vstack((object_info , index))

plot_y = [float(row[1]) for row in object_info]
del plot_y[0]


#Conversions to ABMAG
def ABMAG_Convert(FluxVals):
    ABMAG_list = []
    for i in range(0, len(FluxVals)):
        temp = abs(FluxVals[i]) /4
        ABMAG =  ((-2.5 * np.log(temp)) +23.9)
        ABMAG_list.append(ABMAG)


#Plotting just a couple objects
test_x = []

test_x =  hdulist[0].data[15338:15438 , 8683:8783]
# print(np.shape(test_x))


print(x_coords[9])
print(y_coords[9])
#plt.plot(test_x)
plt.imshow(test_x, cmap = 'gray')
plt.colorbar()
#plt.scatter(plot_y, objects_Imag)
# plt.plot(objects_Imag)
plt.show()



