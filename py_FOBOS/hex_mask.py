def hex_mask(image, x0, y0, r_aper, arr_size, plot = False):
"""
Code that will generate all the pixels of each fiber based on the central pixel, radius, and array size.
Important: The code divides the pixels unevenly and so some fibers have more pixels than other fibers
Note: The outer rings tend to have more pixels than the inner ones
Goal: Preserve as much flux and to not overlap pixels
Important Note: This code does not output the flux of each pixel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
INPUT:
------
     
     image: data.type array usually
         the input image of whatever picture/data image you draw from (Ex: ACS Image)
     x0: int
         x coordinate of the object
     y0: int
         y coordinate of the object
     r_aper: int
         radius of fiber already in units of pixels
     arr_size: int
         the number of Mini-IFU fibers for each object
     plot: true or false
         if you would like to see the plot of the object with the pixels on top of it 

OUTPUT:
-------
     
     pixel_list: dict
         Has subdict x,y which has its own subdicts for each string number of all the pixels
     x_cen: dict
         a dictionary for x coordinates with subdictionaries of each individual fiber_number
     y_cen: dict 
         a dictionary for y coordinates with subdictionaries of each individual fiber_number
     plot(True): 
         plot the image with the pixels of each fiber in color to denote the different fiber
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# ASSUMPTION: The objects of the image will be far enough away from each other 
# so there will be no overlapping of data between each object's fiber
    
from skimage.draw import polygon
from skimage.draw import polygon_perimeter
import numpy as np

"""
---------
    Running fiber_center function that will find the center of each fiber 
    which outputs dictionary x/y_cen[str(fiber_number)] and the bin mask.
    Also generating a dictionary called pixel_list that will be filled subdicts of ['x'] and ['y'] 
    These subdicts will contain dictionaries of each fiber [str(fiber_number)]
    Then each of those wil have all the points
---------
"""
    x_cen, y_cen, img = fiber_centers(x0, y0, r_aper, arr_size)
    pixel_list = {}
    
    #If only looking for one fiber this will output in a slightly different format. 
    if arr_size == 1:
        
        #Calculating the vertices and rounding down so there will not be holes in the pixels later on.
        x = x_cen['1']
        y = y_cen['1']
        height = np.floor(r_aper * np.cos( np.pi/6))
        width = np.floor(r_aper * np.sin( np.pi/6))
        vertex_x, vertex_y = x + width, y + height
        
        #Defining the vertices here and then inputting it into skimage's polygon and polygon perimeter.
        c = np.array([ x - r_aper , (2 * x) - vertex_x , vertex_x , x + r_aper, vertex_x , (2 * x) - vertex_x, x - r_aper])
        r = np.array([ y , vertex_y, vertex_y, y, (2 * y) - vertex_y , (2 * y) - vertex_y , y])
        yy, xx = polygon_perimeter(r,c)
        rr, cc = polygon(yy, xx)
        
        #Saving all the points into a dictionary pixel_list['x/y'][str(fiber_number)]
        rr = np.append([rr],[yy])
        cc = np.append([cc],[xx])
        pixel_list['x'] = cc
        pixel_list['y'] = rr
        img[rr, cc] = 1
    
    #For cases in which number of fibers is greater than 1
    if arr_size > 1:
        pixel_list['x'] = {}
        pixel_list['y'] = {}
        for i in range(0, arr_size):
            
            #Calculating the vertices and rounding down so there will not be holes in the pixels later on.
            x = x_cen[str(i + 1)]
            y = y_cen[str(i + 1)]
            height = np.floor(r_aper * np.cos( np.pi/6))
            width = np.floor(r_aper * np.sin( np.pi/6))
            vertex_x, vertex_y = (x + width), (y + height)
           
            #Defining the vertices here and then inputting it into skimage's polygon and polygon perimeter.
            c = np.array([ x - r_aper, (2 * x) - vertex_x , vertex_x, x + r_aper , vertex_x, (2 * x) - vertex_x , x - r_aper])
            r = np.array([ y , vertex_y, vertex_y, y, (2 * y) - vertex_y , (2 * y) - vertex_y , y])
            yy, xx = polygon_perimeter(r, c)
            rr, cc = polygon(yy, xx)
            
            #Saving all the points into a dictionary pixel_list['x/y'][str(fiber_number)]
            #For the bin mask we add 1 to the previous mask so we can find the overlap later on. 
            rr = np.append([rr],[yy])
            cc = np.append([cc],[xx])
            pixel_list['x'][str(i + 1)] = cc
            pixel_list['y'][str(i + 1)] = rr
            img[rr,cc] = img[rr,cc] + 1
        
        #Now we want to find the overlap in which the bin mask bits are greater than 1
        overlap = np.where(img > 1)
        for i in range(0, len(overlap[0])):
            overlap_difference = []
            
            #Calculating where the nearest fiber is based on the distance from the overlapping tuple to the center of the other fibers.
            for j in range(0, arr_size):
                difference = ((overlap[0][i] - y_cen[str(j + 1)]) ** 2) + ((overlap[1][i] - x_cen[str(j + 1)]) ** 2)
                overlap_difference.append(difference)
            
            #Trying to find the minimum difference (where the objects are closest).
            #Factoring in pixels that might be equidistant to other objects.
            #Assigning pixels to the fiber with fewer pixels and then assigning to the smaller fiber_number
            #We want to have more flux in the central fibers and this should account for that in those cases.
            assignment = [e for e,x in enumerate(overlap_difference) if x == np.min(overlap_difference)]
            if len(assignment) == 1:
                assignment = assignment[0] + 1
            elif len(assignment) > 1:
                if len(assignment) > 2:
                    print(assignment)
                if len(pixel_list['x'][str(assignment[0] + 1)]) > len(pixel_list['x'][str(assignment[1] + 1)]):
                        assignment = assignment[1] + 1
                elif len(pixel_list['x'][str(assignment[0] + 1)]) < len(pixel_list['x'][str(assignment[1] + 1)]):
                        assignment = assignment[0] + 1
                elif len(pixel_list['x'][str(assignment[0] + 1)]) == len(pixel_list['x'][str(assignment[1] + 1)]):
                        if assignment[0] < assignment[1]:
                            assignment = assignment[0] + 1
                        else:
                            assignment = assignment[1] + 1
            
            #Going into each fiber to delete the specific x,y tuple to remove the potential overlap 
            #Then adding the pixel back into the correct fiber number
            for k in range(0, arr_size):
                delete = np.where(np.logical_and(pixel_list['y'][str(k + 1)] == overlap[0][i], pixel_list['x'][str(k + 1)] == overlap[1][i]))
                if len(delete[0]) > 0:
                    pixel_list['y'][str(k+1)] = np.delete(pixel_list['y'][str(k+1)], delete[0])
                    pixel_list['x'][str(k+1)] = np.delete(pixel_list['x'][str(k+1)], delete[0])
            pixel_list['y'][str(assignment)] = np.append(pixel_list['y'][str(assignment)], overlap[0][i])
            pixel_list['x'][str(assignment)] = np.append(pixel_list['x'][str(assignment)], overlap[1][i])
            img[pixel_list['y'][str(assignment)][-1],pixel_list['x'][str(assignment)][-1]] = 1
    
    #Plot if true and return the pixels
    if plot == False:
        return pixel_list, x_cen, y_cen, 
    if plot == True:
        for i in range(0, arr_size):
            plt.scatter(pixel_list['y'][str(i + 1)], pixel_list['x'][str(i + 1)])
        plt.imshow(image)
        plt.show()
        return pixel_list, x_cen, y_cen

def fiber_centers(image, x0, y0, r_aper, arr_size):
"""
Generating fiber centers by taking the central pixel and using basic trigonometry to calculate vertices
Important: Since there might be weird decimals the code rounds down
Also: There must be symmetry which might be lost in the decimals 
so the code reflects each point by calculating differences

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
INPUT:
------
     
     image: data.type array usually
         the input image of whatever picture/data image you draw from (Ex: ACS Image)
     x0: int
         x coordinate of the object
     y0: int
         y coordinate of the object
     r_aper: int
         radius of fiber already in units of pixels
     arr_size: int
         the number of Mini-IFU fibers for each object


OUTPUT:
-------
     
     pixel_list: dict
         Has subdict x,y which has its own subdicts for each string number of all the pixels
     x_cen: dict
         a dictionary for x coordinates with subdictionaries of each individual fiber_number
     y_cen: dict 
         a dictionary for y coordinates with subdictionaries of each individual fiber_number
     img: array
         bin mask of 1's and 0's where 1 is a positive hit on a x and y coordinate of a pixel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

        
        #Generating the dicts and setting up the bin mask
        #Bin mask will be the size of shape of the image
        x_cen = {}
        y_cen = {}
        x , y = np.shape(image)
        img = np.zeros((y, x), dtype=np.uint16)
        
        #Adding the central pixel to the dict and the format will stay consistent throughout
        #Finding the center of each fiber using some basic trig
        x_cen['1'], y_cen['1'] = x0, y0
        if arr_size > 1:
            x_cen['2'], y_cen['2'] = (x0 - r_aper - ((r_aper * np.sin(np.pi/6)))), (y0 + np.floor(r_aper * np.cos(np.pi / 6)))
            x_cen['3'], y_cen['3'] =  x0, (y0 - 1 + np.floor(2 * (r_aper * np.cos(np.pi/6))))
            x_cen['4'], y_cen['4'] = (2* x0) - x_cen['2'], y_cen['2']
            x_cen['5'], y_cen['5'] = x_cen['4'], (2 * y0) - y_cen['4']
            x_cen['6'], y_cen['6'] = x_cen['3'], (2 * y0) - y_cen['3']
            x_cen['7'], y_cen['7'] = x_cen['2'], (2 * y0) - y_cen['2']
            if arr_size > 7:
                x_cen['8'], y_cen['8'] = x0 - (2 * (x_cen['2'] - x0)), y0
                x_cen['9'], y_cen['9'] = x_cen['8'], y_cen['3']
                x_cen['10'], y_cen['10'] = x_cen['2'], (2 * y_cen['3']) - y_cen['2']
                x_cen['11'], y_cen['11'] = x0, y0 +( 2 * (y_cen['3'] - y0))
                x_cen['12'], y_cen['12'] = x_cen['4'], y_cen['10']
                x_cen['13'], y_cen['13'] = x0 + (x0 - x_cen['8']), y_cen['3']
                x_cen['14'], y_cen['14'] = x_cen['13'], y0
                x_cen['15'], y_cen['15'] = x_cen['13'], y_cen['6']
                x_cen['16'], y_cen['16'] = x_cen['4'], (2 * y0) - y_cen['12'] 
                x_cen['17'], y_cen['17'] = x0, y0 - (y_cen['11'] - y0)
                x_cen['18'], y_cen['18'] = x_cen['2'], y_cen['16']
                x_cen['19'], y_cen['19'] = x_cen['8'], y_cen['6']
        
        #In the case of one fiber there is only one center so bin mask only has one hit
        #Otherwise each specific fiber has a 1 for the center and we've populated the bin mask. 
        if arr_size == 1:
            img[y_cen['1'], x_cen['1']] = 1
        else:
            for i in range(0, arr_size):
                    x_cen[str(i + 1)] = np.int(x_cen[str(i + 1)])
                    y_cen[str(i + 1)] = np.int(y_cen[str(i + 1)])
                    img[y_cen[str(i + 1)], x_cen[str(i + 1)]] = 1
        plt.imshow(img)
        plt.show()
        return x_cen, y_cen, img
