import numpy as np
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import scipy.ndimage.filters as ndimage
class gen_spectra:
        def __init__(self, fits_tile):
                self.hdulist = fits.open(fits_tile)
                #In degrees!
                RA = self.hdulist[0].header['RA_Targ']
                DEC = self.hdulist[0].header['DEC_Targ']

                #Arcsecond:Pixel Conversion factor and Flux Area
                #.57 arcseconds is .015833333 degrees
                self.Dia = .57
                self.radi_arcs = self.Dia/2
                self.outer_radi =.01583333333 * np.sqrt(7)
                self.A_P = .03
                self.radi = ((self.Dia/self.A_P)/2)
                xy_max = int(self.radi/(np.sqrt(2)))
                self.radi_outer = self.radi * np.sqrt(7)

                #Target_index will be 0 for now in case we aren't looking at specific targets. 
                self.target_index = 0
                #Each pixel
                self.scidata = self.hdulist[0].data

                #Instead convert the 4 points of the box to create the dimensions
                
                self.box_4= wcs.WCS(self.hdulist[0].header)
                x_range, y_range = np.shape(self.scidata)
                box_coords = np.array([(0 + xy_max, 0  + xy_max),(0 + xy_max, y_range - xy_max),(x_range - xy_max, 0 + xy_max),(x_range - xy_max, y_range - xy_max)])
                Box_RA_DEC = self.box_4.wcs_pix2world(box_coords, 1)
                RA_m = Box_RA_DEC[2,0].tolist()
                RA_ma = Box_RA_DEC[0,0].tolist()
                DEC_m = Box_RA_DEC[0,1].tolist()
                DEC_ma = Box_RA_DEC[1,1].tolist()
                self.RA_min = RA_m + self.outer_radi
                self.RA_max = RA_ma - self.outer_radi
                self.DEC_min = DEC_m + self.outer_radi
                self.DEC_max = DEC_ma - self.outer_radi


#From the catalog find all the targets within the range of the image
#Generate a list of of all the RA and DEC
#Using if statement look through each ID to find values in the range
#of the RA and DEC of the specific .fits file
#Finding only the relevant catalog_RA and catalog_DEC

        def catalog_info(self, fits_catalog):
                self.catalog = fits.open(fits_catalog)
                self.catalog_wcs = wcs.WCS(self.catalog[1].header)
                self.catalog_id = self.catalog[1].data['id']
                self.catalog_RA = np.array(self.catalog[1].data['RA'])
                self.catalog_DEC = np.array(self.catalog[1].data['DEC'])
                self.catalog_Imag = np.array(self.catalog[1].data['imag'])
                self.catalog_Umag = np.array(self.catalog[1].data['umag'])
                self.catalog_Zmag = np.array(self.catalog[1].data['zmag'])
                self.catalog_Rmag = np.array(self.catalog[1].data['rmag'])
                self.catalog_Vmag = np.array(self.catalog[1].data['vmag'])
                self.catalog_Gmag = np.array(self.catalog[1].data['gmag'])
                self.catalog_Bmag = np.array(self.catalog[1].data['bmag'])
                self.catalog_zgal = np.array(self.catalog[1].data['zgal'])
                self.catalog_zilb = np.array(self.catalog[1].data['z_ilb'])
                self.catalog_chis = np.array(self.catalog[1].data['chi_sec'])
                self.catalog_zsec = np.array(self.catalog[1].data['zsec'])

                self.objects_index = np.where(np.logical_and(self.catalog_RA  <= self.RA_max , self.catalog_RA >= self.RA_min))
                self.objects_index = np.intersect1d(self.objects_index, np.where(np.logical_and(self.catalog_DEC <= self.DEC_max , self.catalog_DEC >= self.DEC_min)))
                self.objects_index = np.intersect1d(self.objects_index, np.where(np.logical_and(self.catalog_Imag > 0 ,self.catalog_Umag > 0)))
                self.objects_index = np.intersect1d(self.objects_index, np.where(np.logical_and(self.catalog_Zmag > 0 ,self.catalog_Vmag > 0)))
                self.objects_index = np.intersect1d(self.objects_index, np.where(np.logical_and(self.catalog_Rmag > 0 ,self.catalog_Gmag > 0)))
                self.objects_index = np.intersect1d(self.objects_index, np.where(np.logical_and(self.catalog_zgal > 0 ,self.catalog_zilb > 0)))
                self.objects_index = np.intersect1d(self.objects_index, np.where(np.logical_and(self.catalog_Bmag > 0 ,self.catalog_chis > 0)))
                self.objects_id = self.catalog_id[self.objects_index]
                self.objects_RA = self.catalog_RA[self.objects_index]
                self.objects_DEC = self.catalog_DEC[self.objects_index]
                self.objects_Imag = self.catalog_Imag[self.objects_index]
                self.objects_Umag = self.catalog_Umag[self.objects_index]
                self.objects_Zmag = self.catalog_Zmag[self.objects_index]
                self.objects_Vmag = self.catalog_Vmag[self.objects_index]
                self.objects_Rmag = self.catalog_Rmag[self.objects_index]
                self.objects_Gmag = self.catalog_Gmag[self.objects_index]
                self.objects_Bmag =  self.catalog_Bmag[self.objects_index]


                #Most x,y coords were in decimals so I made them into integers. Some of them rounded upwards
                x = np.array([(self.objects_RA[0], self.objects_DEC[0])], np.float)
                for i in range(1, len(self.objects_index)):
                        y = np.array([(self.objects_RA[i], self.objects_DEC[i])])
                        x = np.vstack(( x , y))
                self.objects = (self.box_4.wcs_world2pix(x , 1))
                self.objects = self.objects.astype(int)
                self.x_coords = [row[0] for row in self.objects]
                self.y_coords = [row[1] for row in self.objects]
        def target_image(self, fits_plot):
                self.box_image = fits_plot[(self.y_coords[self.target_index] - 50):(self.y_coords[self.target_index] + 50), (self.x_coords[self.target_index] - 50): (self.x_coords[self.target_index] + 50)]
        def overplotting(self):
                #Creating the first inner fiber
                circle_row = []
                circle_column = []

                for i in range(1 , 10):
                        for j in range(1 , 10):
                                if 0 <= (i ** 2) + (j ** 2) <= (self.radi ** 2):
                                        circle_row.append(i)
                                        circle_column.append(j)
                for i in range(len(circle_row)):
                    circle_row.append( circle_row[i])
                    circle_column.append(-circle_column[i])
                for j in range(len(circle_row)):
                        circle_row.append(-circle_row[j])
                        circle_column.append(circle_column[j])

                for i in range(1,10):
                        if 0 <= (i ** 2) <= (self.radi ** 2):
                                circle_row.append(i)
                                circle_column.append(0)
                                circle_row.append(0)
                                circle_column.append(i)
                                circle_row.append(-i)
                                circle_column.append(0)
                                circle_row.append(0)
                                circle_column.append(-i)

                #Generating all the x,y coordinates for each object
                self.object_x_coords = []
                self.object_y_coords = []
                self.fiber0_xcoords = []
                self.fiber0_ycoords = []

                for object_number in range(0, len(self.x_coords)):
                        fiber_row = [self.y_coords[object_number] + y for y in circle_row]
                        fiber_column = [self.x_coords[object_number] + x for x in circle_column]
                        self.fiber0_xcoords.append(self.x_coords[object_number])
                        self.fiber0_ycoords.append(self.y_coords[object_number])
                        self.fiber0_xcoords.extend(fiber_column)
                        self.fiber0_ycoords.extend(fiber_row)

                #NOTE: EACH OBJECT has 293 pixels for fiber0
                #Generate flux values for each object's fiber0

                self.flux0_values = []
                for i in range(0,len(self.fiber0_xcoords)):
                        self.flux0_values.append(self.scidata[self.fiber0_ycoords[i] , self.fiber0_xcoords[i]])
                #Radius of the outer circle
                #All the points on the outer circle

                outer_circle_row = []
                outer_circle_column = []
                for i in range(1 , 26):
                        for j in range(1 , 26):
                                if (self.radi ** 2)  < (i ** 2) + (j ** 2) <= (self.radi_outer ** 2):
                                        outer_circle_row.append(i)
                                        outer_circle_column.append(j)

                for i in range(len(outer_circle_row)):
                    outer_circle_row.append( outer_circle_row[i])
                    outer_circle_column.append(-outer_circle_column[i])
                for j in range(len(outer_circle_row)):
                        outer_circle_row.append(-outer_circle_row[j])
                        outer_circle_column.append(outer_circle_column[j])

                for i in range(0,26):
                        if (self.radi ** 2) < (i ** 2) <= (self.radi_outer ** 2):
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

                self.fiber1_x = []
                self.fiber1_y = []
                self.fiber2_x = []
                self.fiber2_y = []
                self.fiber3_x = []
                self.fiber3_y = []
                self.fiber4_x = []
                self.fiber4_y = []
                self.fiber5_x = []
                self.fiber5_y = []
                self.fiber6_x = []
                self.fiber6_y = []

                for i in range(0 , len(outer_circle_row)):
                        if outer_circle_row[i] > 0 and outer_circle_column[i] > 0:
                                if 0  <  (outer_circle_row[i]/outer_circle_column[i])  < 1.65:
                                        self.fiber1_x.append(outer_circle_column[i])
                                        self.fiber1_y.append(outer_circle_row[i])
                                else:
                                        self.fiber2_x.append(outer_circle_column[i])
                                        self.fiber2_y.append(outer_circle_row[i])
                        if outer_circle_row[i] < 0 and outer_circle_column[i] > 0:
                                if 0  <  (-outer_circle_row[i]/outer_circle_column[i])  <= 1.85:
                                        self.fiber6_x.append(outer_circle_column[i])
                                        self.fiber6_y.append(outer_circle_row[i])
                                else:
                                        self.fiber5_x.append(outer_circle_column[i])
                                        self.fiber5_y.append(outer_circle_row[i])
                        if outer_circle_row[i] > 0 and outer_circle_column[i] < 0:
                                if 0  <  (outer_circle_row[i]/-outer_circle_column[i])  <= 1.85:
                                        self.fiber3_x.append(outer_circle_column[i])
                                        self.fiber3_y.append(outer_circle_row[i])
                                else:
                                        self.fiber2_x.append(outer_circle_column[i])
                                        self.fiber2_y.append(outer_circle_row[i])
                        if outer_circle_row[i] < 0 and outer_circle_column[i] < 0:
                                if 0  <  (-outer_circle_row[i]/-outer_circle_column[i])  < 1.65:
                                        self.fiber4_x.append(outer_circle_column[i])
                                        self.fiber4_y.append(outer_circle_row[i])
                                else:
                                        self.fiber5_x.append(outer_circle_column[i])
                                        self.fiber5_y.append(outer_circle_row[i])
                for i in range(0 , len(outer_circle_row)):
                        if outer_circle_row[i] == 0 and outer_circle_column[i] > 0:
                                self.fiber1_x.append(outer_circle_column[i])
                                self.fiber1_y.append(outer_circle_row[i])
                        if outer_circle_row[i] > 0 and outer_circle_column[i] == 0:
                                self.fiber2_x.append(outer_circle_column[i])
                                self.fiber2_y.append(outer_circle_row[i])
                        if outer_circle_row[i] == 0 and outer_circle_column[i] < 0:
                                self.fiber4_x.append(outer_circle_column[i])
                                self.fiber4_y.append(outer_circle_row[i])
                        if outer_circle_row[i] < 0 and outer_circle_column[i] == 0:
                                self.fiber5_x.append(outer_circle_column[i])
                                self.fiber5_y.append(outer_circle_row[i])

                #Creating fibers for each object
                self.fiber1_xcoords = []
                self.fiber1_ycoords = []
                self.fiber2_xcoords = []
                self.fiber2_ycoords = []
                self.fiber3_xcoords = []
                self.fiber3_ycoords = []
                self.fiber4_xcoords = []
                self.fiber4_ycoords = []
                self.fiber5_xcoords = []
                self.fiber5_ycoords = []
                self.fiber6_xcoords = []
                self.fiber6_ycoords = []


                for object_number in range(0 , len(self.y_coords)):
                        outer_fiber1_row = [self.y_coords[object_number] + y for y in self.fiber1_y]
                        outer_fiber1_column = [self.x_coords[object_number] + x for x in self.fiber1_x]
                        outer_fiber2_row = [self.y_coords[object_number] + y for y in self.fiber2_y]
                        outer_fiber2_column = [self.x_coords[object_number] + x for x in self.fiber2_x]
                        outer_fiber3_row = [self.y_coords[object_number] + y for y in self.fiber3_y]
                        outer_fiber3_column = [self.x_coords[object_number] + x for x in self.fiber3_x]
                        outer_fiber4_row = [self.y_coords[object_number] + y for y in self.fiber4_y]
                        outer_fiber4_column = [self.x_coords[object_number] + x for x in self.fiber4_x]
                        outer_fiber5_row = [self.y_coords[object_number] + y for y in self.fiber5_y]
                        outer_fiber5_column = [self.x_coords[object_number] + x for x in self.fiber5_x]
                        outer_fiber6_row = [self.y_coords[object_number] + y for y in self.fiber6_y]
                        outer_fiber6_column = [self.x_coords[object_number] + x for x in self.fiber6_x]   
                        self.fiber1_xcoords.extend(outer_fiber1_column)
                        self.fiber1_ycoords.extend(outer_fiber1_row)
                        self.fiber2_xcoords.extend(outer_fiber2_column)
                        self.fiber2_ycoords.extend(outer_fiber2_row)
                        self.fiber3_xcoords.extend(outer_fiber3_column)
                        self.fiber3_ycoords.extend(outer_fiber3_row)
                        self.fiber4_xcoords.extend(outer_fiber4_column)
                        self.fiber4_ycoords.extend(outer_fiber4_row)
                        self.fiber5_xcoords.extend(outer_fiber5_column)
                        self.fiber5_ycoords.extend(outer_fiber5_row)
                        self.fiber6_xcoords.extend(outer_fiber6_column)
                        self.fiber6_ycoords.extend(outer_fiber6_row)


                #Given pixel coordinates I set new variables for each individual pixel's flux value
                #Then I will create a new list of all the objects with variables set to their self.fiber's total flux
                #Fiber1/4: 285
                #Fiber2/5: 282
                #Fiber3/6: 283

                self.flux1_values = []
                self.flux2_values = []
                self.flux3_values = []
                self.flux4_values = []
                self.flux5_values = []
                self.flux6_values = []

                for i in range(0,len(self.fiber1_xcoords)):
                        self.flux1_values.append(self.scidata[self.fiber1_ycoords[i] , self.fiber1_xcoords[i]])
                for i in range(0,len(self.fiber2_xcoords)):
                        self.flux2_values.append(self.scidata[self.fiber2_ycoords[i] , self.fiber2_xcoords[i]])
                for i in range(0,len(self.fiber3_xcoords)):
                        self.flux3_values.append(self.scidata[self.fiber3_ycoords[i] , self.fiber3_xcoords[i]])
                for i in range(0,len(self.fiber4_xcoords)):
                        self.flux4_values.append(self.scidata[self.fiber4_ycoords[i] , self.fiber4_xcoords[i]])
                for i in range(0,len(self.fiber5_xcoords)):
                        self.flux5_values.append(self.scidata[self.fiber5_ycoords[i] , self.fiber5_xcoords[i]])
                for i in range(0,len(self.fiber6_xcoords)):
                        self.flux6_values.append(self.scidata[self.fiber6_ycoords[i] , self.fiber6_xcoords[i]])
#Summing all the points for each Mini-IFU
        def Flux_Sum(self):
                self.IFU_flux0 = []
                self.IFU_flux1 = []
                self.IFU_flux2 = []
                self.IFU_flux3 = []
                self.IFU_flux4 = []
                self.IFU_flux5 = []
                self.IFU_flux6 = []
                self.fiber_sum = {}
                for i in range(0, len(self.objects_id)):
                        counter0 = i * 293
                        counter1 = i * 285
                        counter2 = i * 282
                        counter3 = i * 283
                        self.fiber_sum[str(self.objects_id[i]) + '_fiber0'.format(i)] = sum(self.flux0_values[counter0 + j] for j in range(293))
                        self.fiber_sum[str(self.objects_id[i]) + '_fiber1'.format(i)] = sum(self.flux1_values[counter1 + k] for k in range(285))
                        self.fiber_sum[str(self.objects_id[i]) + '_fiber2'.format(i)] = sum(self.flux2_values[counter2 + l] for l in range(282))
                        self.fiber_sum[str(self.objects_id[i]) + '_fiber3'.format(i)] = sum(self.flux3_values[counter3 + m] for m in range(283))
                        self.fiber_sum[str(self.objects_id[i]) + '_fiber4'.format(i)] = sum(self.flux4_values[counter1 + n] for n in range(285))
                        self.fiber_sum[str(self.objects_id[i]) + '_fiber5'.format(i)] = sum(self.flux5_values[counter2 + o] for o in range(282))
                        self.fiber_sum[str(self.objects_id[i]) + '_fiber6'.format(i)] = sum(self.flux6_values[counter3 + p] for p in range(283))
                        self.IFU_flux0.append(self.fiber_sum[str(self.objects_id[i]) + '_fiber0'])
                        self.IFU_flux1.append(self.fiber_sum[str(self.objects_id[i]) + '_fiber1'])
                        self.IFU_flux2.append(self.fiber_sum[str(self.objects_id[i]) + '_fiber2'])
                        self.IFU_flux3.append(self.fiber_sum[str(self.objects_id[i]) + '_fiber3'])
                        self.IFU_flux4.append(self.fiber_sum[str(self.objects_id[i]) + '_fiber4'])
                        self.IFU_flux5.append(self.fiber_sum[str(self.objects_id[i]) + '_fiber5'])
                        self.IFU_flux6.append(self.fiber_sum[str(self.objects_id[i]) + '_fiber6'])
#Integrate Flux values
        def int_flux(self):
                self.IFU_sum = []
                for i in range(0, len(self.objects_id)):
                        flux_sum = [float(self.fiber_sum[str(self.objects_id[i]) + '_fiber0']) , float(self.fiber_sum[str(self.objects_id[i]) + '_fiber1']), float(self.fiber_sum[str(self.objects_id[i]) + '_fiber2']) , float(self.fiber_sum[str(self.objects_id[i]) + '_fiber3']) , float(self.fiber_sum[str(self.objects_id[i]) + '_fiber4']) , float(self.fiber_sum[str(self.objects_id[i]) + '_fiber5']) , float(self.fiber_sum[str(self.objects_id[i]) + '_fiber6'])]
                        self.IFU_sum.append(sum(flux_sum))

#Plotting just a couple objects
        def graph_image(self, image):
                plt.imshow(image, cmap = 'gray')
                plt.gca().invert_yaxis()
                plt.colorbar()
                plt.show()
        def ABMAG_Convert(self, FluxVals):
                if isinstance(FluxVals, str):
                        FluxVals = eval(FluxVals)
                self.ABMAG_list = []
                for i in range(0, len(FluxVals)):
                        ABMAG =  ((-2.5 * np.log10(FluxVals[i])) +23.9)
                        self.ABMAG_list.append(ABMAG)
        def Flux_Convert(self, ABMAGVals):
                if isinstance(ABMAGVals, str):
                        ABMAGVals = eval(ABMAGVals)
                self.Flux_list = []
                for i in range(0, len(ABMAGVals)):
                        FluxVals = 10 ** ((ABMAGVals[i] - 23.9)/(-2.5))
                        self.Flux_list.append(FluxVals)
        def border_plot(self):
                #Generating the border for the inner Mini-IFU
                border_row = []
                border_column = []

                for i in range(1 , 10):
                        for j in range(1 , 10):
                                if 75 <= (i ** 2) + (j ** 2) <= (self.radi ** 2):
                                        border_row.append(i)
                                        border_column.append(j)
                for i in range(len(border_row)):
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
                                if 600 < (i ** 2) + (j ** 2) <= (self.radi_outer ** 2):
                                        outer_border_row.append(i)
                                        outer_border_column.append(j)
                for i in range(9 , 23):
                        for j in range(1 , 23):
                                if -.5 < i - (j* np.sqrt(3)) < .5:
                                        outer_border_row.append(i)
                                        outer_border_column.append(j)


                for i in range(len(outer_border_row)):
                    outer_border_row.append( outer_border_row[i])
                    outer_border_column.append(-outer_border_column[i])
                for j in range(len(outer_border_row)):
                        outer_border_row.append(-outer_border_row[j])
                        outer_border_column.append(outer_border_column[j])

                for i in range(0,26):
                        if (self.radi ** 2) < (i ** 2) <= (self.radi_outer ** 2):
                                outer_border_row.append(0)
                                outer_border_column.append(i)
                                outer_border_row.append(0)
                                outer_border_column.append(-i)


                #Dividing the outer border to respective fibers
                self.border_fiber1_x = []
                self.border_fiber1_y = []
                self.border_fiber2_x = []
                self.border_fiber2_y = []
                self.border_fiber3_x = []
                self.border_fiber3_y = []
                self.border_fiber4_x = []
                self.border_fiber4_y = []
                self.border_fiber5_x = []
                self.border_fiber5_y = []
                self.border_fiber6_x = []
                self.border_fiber6_y = []

                for i in range(0 , len(outer_border_row)):
                        if outer_border_row[i] > 0 and outer_border_column[i] > 0:
                                if 0  <  (outer_border_row[i]/outer_border_column[i])  < 1.65:
                                        self.border_fiber1_x.append(outer_border_column[i])
                                        self.border_fiber1_y.append(outer_border_row[i])
                                else:
                                        self.border_fiber2_x.append(outer_border_column[i])
                                        self.border_fiber2_y.append(outer_border_row[i])
                        if outer_border_row[i] < 0 and outer_border_column[i] > 0:
                                if 0  <  (-outer_border_row[i]/outer_border_column[i])  <= 1.85:
                                        self.border_fiber6_x.append(outer_border_column[i])
                                        self.border_fiber6_y.append(outer_border_row[i])
                                else:
                                        self.border_fiber5_x.append(outer_border_column[i])
                                        self.border_fiber5_y.append(outer_border_row[i])
                        if outer_border_row[i] > 0 and outer_border_column[i] < 0:
                                if 0  <  (outer_border_row[i]/-outer_border_column[i])  <= 1.85:
                                        self.border_fiber3_x.append(outer_border_column[i])
                                        self.border_fiber3_y.append(outer_border_row[i])
                                else:
                                        self.border_fiber2_x.append(outer_border_column[i])
                                        self.border_fiber2_y.append(outer_border_row[i])
                        if outer_border_row[i] < 0 and outer_border_column[i] < 0:
                                if 0  <  (-outer_border_row[i]/-outer_border_column[i])  < 1.65:
                                        self.border_fiber4_x.append(outer_border_column[i])
                                        self.border_fiber4_y.append(outer_border_row[i])
                                else:
                                        self.border_fiber5_x.append(outer_border_column[i])
                                        self.border_fiber5_y.append(outer_border_row[i])
                for i in range(0 , len(outer_border_row)):
                        if outer_border_row[i] == 0 and outer_border_column[i] > 0:
                                self.border_fiber1_x.append(outer_border_column[i])
                                self.border_fiber1_y.append(outer_border_row[i])
                        if outer_border_row[i] > 0 and outer_border_column[i] == 0:
                                self.border_fiber2_x.append(outer_border_column[i])
                                self.border_fiber2_y.append(outer_border_row[i])
                        if outer_border_row[i] == 0 and outer_border_column[i] < 0:
                                self.border_fiber4_x.append(outer_border_column[i])
                                self.border_fiber4_y.append(outer_border_row[i])
                        if outer_border_row[i] < 0 and outer_border_column[i] == 0:
                                self.border_fiber5_x.append(outer_border_column[i])
                                self.border_fiber5_y.append(outer_border_row[i])

                #Creating fibers for each object
                self.border_fiber1_xcoords = []
                self.border_fiber1_ycoords = []
                self.border_fiber2_xcoords = []
                self.border_fiber2_ycoords = []
                self.border_fiber3_xcoords = []
                self.border_fiber3_ycoords = []
                self.border_fiber4_xcoords = []
                self.border_fiber4_ycoords = []
                self.border_fiber5_xcoords = []
                self.border_fiber5_ycoords = []
                self.border_fiber6_xcoords = []
                self.border_fiber6_ycoords = []
                border_fiber1_row = [50 + y for y in self.border_fiber1_y]
                border_fiber1_column = [50 + x for x in self.border_fiber1_x]
                border_fiber2_row = [50 + y for y in self.border_fiber2_y]
                border_fiber2_column = [50 + x for x in self.border_fiber2_x]
                border_fiber3_row = [50 + y for y in self.border_fiber3_y]
                border_fiber3_column = [50 + x for x in self.border_fiber3_x]
                border_fiber4_row = [50 + y for y in self.border_fiber4_y]
                border_fiber4_column = [50 + x for x in self.border_fiber4_x]
                border_fiber5_row = [50 + y for y in self.border_fiber5_y]
                border_fiber5_column = [50 + x for x in self.border_fiber5_x]
                border_fiber6_row = [50 + y for y in self.border_fiber6_y]
                border_fiber6_column = [50 + x for x in self.border_fiber6_x]  
                self.border_fiber1_xcoords.extend(border_fiber1_column)
                self.border_fiber1_ycoords.extend(border_fiber1_row)
                self.border_fiber2_xcoords.extend(border_fiber2_column)
                self.border_fiber2_ycoords.extend(border_fiber2_row)
                self.border_fiber3_xcoords.extend(border_fiber3_column)
                self.border_fiber3_ycoords.extend(border_fiber3_row)
                self.border_fiber4_xcoords.extend(border_fiber4_column)
                self.border_fiber4_ycoords.extend(border_fiber4_row)
                self.border_fiber5_xcoords.extend(border_fiber5_column)
                self.border_fiber5_ycoords.extend(border_fiber5_row)
                self.border_fiber6_xcoords.extend(border_fiber6_column)
                self.border_fiber6_ycoords.extend(border_fiber6_row)

                #Generating all the x,y coordinates for each object

                self.border_fiber0_xcoords = []
                self.border_fiber0_ycoords = []
                fiber0_border_row = [ y + 50  for y in border_row]
                fiber0_border_column = [ x + 50 for x in border_column]
                self.border_fiber0_xcoords.append(self.x_coords[self.target_index])
                self.border_fiber0_ycoords.append(self.y_coords[self.target_index])
                self.border_fiber0_xcoords.extend(fiber0_border_column)
                self.border_fiber0_ycoords.extend(fiber0_border_row)
        def fwhm2sigma(self, fwhm):
                self.fwhm = fwhm
                return fwhm / np.sqrt(8 * np.log(2))
        def smooth(self, fwhm):
                self.fwhm = fwhm
                self.sigma1 = self.fwhm2sigma(fwhm)/self.A_P
                self.smoothed_image = ndimage.gaussian_filter(self.scidata, sigma = self.sigma1)
                self.target_image(self.smoothed_image)
        def smoothing(self, fwhm):
                self.smoothed_IFU_flux = []
                self.smoothed_IFU_sum = []
                gaussian_0 = ndimage.gaussian_filter(self.flux0_values, sigma = self.sigma1)
                gaussian_1 = ndimage.gaussian_filter(self.flux1_values, sigma = self.sigma1)
                gaussian_2 = ndimage.gaussian_filter(self.flux2_values, sigma = self.sigma1)
                gaussian_3 = ndimage.gaussian_filter(self.flux3_values, sigma = self.sigma1)
                gaussian_4 = ndimage.gaussian_filter(self.flux4_values, sigma = self.sigma1)
                gaussian_5 = ndimage.gaussian_filter(self.flux5_values, sigma = self.sigma1)
                gaussian_6 = ndimage.gaussian_filter(self.flux6_values, sigma = self.sigma1)

                smooth_sum = {}
                self.smoothed_IFU_flux0 = []
                self.smoothed_IFU_flux1 = []
                self.smoothed_IFU_flux2 = []
                self.smoothed_IFU_flux3 = []
                self.smoothed_IFU_flux4 = []
                self.smoothed_IFU_flux5 = []
                self.smoothed_IFU_flux6 = []
                for i in range(0, len(self.objects_index)):
                        counter0 = i * 293
                        counter1 = i * 285
                        counter2 = i * 282
                        counter3 = i * 283
                        smooth_sum[str(self.objects_id[i]) + '_fiber0'.format(i)] = sum(gaussian_0[counter0 + j] for j in range(293))
                        smooth_sum[str(self.objects_id[i]) + '_fiber1'.format(i)] = sum(gaussian_1[counter1 + k] for k in range(285))
                        smooth_sum[str(self.objects_id[i]) + '_fiber2'.format(i)] = sum(gaussian_2[counter2 + l] for l in range(282))
                        smooth_sum[str(self.objects_id[i]) + '_fiber3'.format(i)] = sum(gaussian_3[counter3 + m] for m in range(283))
                        smooth_sum[str(self.objects_id[i]) + '_fiber4'.format(i)] = sum(gaussian_4[counter1 + n] for n in range(285))
                        smooth_sum[str(self.objects_id[i]) + '_fiber5'.format(i)] = sum(gaussian_5[counter2 + o] for o in range(282))
                        smooth_sum[str(self.objects_id[i]) + '_fiber6'.format(i)] = sum(gaussian_6[counter3 + p] for p in range(283))

                        self.smoothed_IFU_flux0.append(float(smooth_sum[str(self.objects_id[i]) + '_fiber0']))
                        self.smoothed_IFU_flux1.append(float(smooth_sum[str(self.objects_id[i]) + '_fiber1']))
                        self.smoothed_IFU_flux2.append(float(smooth_sum[str(self.objects_id[i]) + '_fiber2']))
                        self.smoothed_IFU_flux3.append(float(smooth_sum[str(self.objects_id[i]) + '_fiber3']))
                        self.smoothed_IFU_flux4.append(float(smooth_sum[str(self.objects_id[i]) + '_fiber4']))
                        self.smoothed_IFU_flux5.append(float(smooth_sum[str(self.objects_id[i]) + '_fiber5']))
                        self.smoothed_IFU_flux6.append(float(smooth_sum[str(self.objects_id[i]) + '_fiber6']))

                        flux_sum = [float(smooth_sum[str(self.objects_id[i]) + '_fiber0']) , float(smooth_sum[str(self.objects_id[i]) + '_fiber1']), float(smooth_sum[str(self.objects_id[i]) + '_fiber2']) , float(smooth_sum[str(self.objects_id[i]) + '_fiber3']) , float(smooth_sum[str(self.objects_id[i]) + '_fiber4']) , float(smooth_sum[str(self.objects_id[i]) + '_fiber5']) , float(smooth_sum[str(self.objects_id[i]) + '_fiber6'])]
                        self.smoothed_IFU_sum.append(float(sum(flux_sum)))
                        
        def overplot_image(self, image):
                plt.scatter(self.border_fiber0_xcoords, self.border_fiber0_ycoords)
                plt.scatter(self.border_fiber1_xcoords, self.border_fiber1_ycoords)
                plt.scatter(self.border_fiber2_xcoords, self.border_fiber2_ycoords)
                plt.scatter(self.border_fiber3_xcoords, self.border_fiber3_ycoords)
                plt.scatter(self.border_fiber4_xcoords, self.border_fiber4_ycoords)
                plt.scatter(self.border_fiber5_xcoords, self.border_fiber5_ycoords)
                plt.scatter(self.border_fiber6_xcoords, self.border_fiber6_ycoords)
                plt.imshow(image , cmap = 'gray')
                plt.gca().invert_yaxis()
                plt.colorbar()
                plt.show()
        # def smooth_flux(self):
        #         print(str(self.fwhm) + "_Smoothed_IFU_Flux0_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_flux0[self.target_index]))
        #         print(str(self.fwhm) + "_Smoothed_IFU_Flux1_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_flux1[self.target_index]))
        #         print(str(self.fwhm) + "_Smoothed_IFU_Flux2_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_flux2[self.target_index]))
        #         print(str(self.fwhm) + "_Smoothed_IFU_Flux3_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_flux3[self.target_index]))
        #         print(str(self.fwhm) + "_Smoothed_IFU_Flux4_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_flux4[self.target_index]))
        #         print(str(self.fwhm) + "_Smoothed_IFU_Flux5_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_flux5[self.target_index]))
        #         print(str(self.fwhm) + "_Smoothed_IFU_Flux6_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_flux6[self.target_index]))
        #         print(str(self.fwhm) + "_Smoothed_IFU_sum_" + str(self.objects_id[self.target_index]) + ": " +  str(self.smoothed_IFU_sum[self.target_index]))
        def normalization(self, mag1, flux1, flux2):
                self.mag2 = mag1 + (-2.5 *np.log10(flux1/flux2))
        def renormalization(self, mag1, mag2, flux2):
                self.flux1val = (flux2 * (10 ** ((mag1 - mag2)/ -2.5)))



               




