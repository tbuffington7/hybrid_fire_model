import yaml
import os
import subprocess
import pandas as pd
import numpy as np
import pdb

class fireLocError(Exception):
    pass


def firecalc(configlocation, room, layout, fireloc):
    #Setting the location of config file, which holds settings
    configs = yaml.load(open(configlocation))
    #In that file, the 'firecalc' section is relevant for this code
    configuse = configs['firecalc']
    #Creating an array of times spaced by timestep that goes to simtime
    timelist = np.arange(0,configuse['simtime'], configuse['timestep'])
    #Importing item info as pandas df
    iteminfo = pd.read_csv(configuse['itemloc'],index_col='item')
    #Importing room specs as pandas df
    roominfo = pd.read_csv(configuse['roomloc'] \
    + '/' + room + '.csv')

    

    #Need room surface area for flashover correlation
    surf_area = 2*roominfo['width'][0]*roominfo['length'][0] \
    + 2*roominfo['width'][0]*roominfo['height'][0] \
    + 2*roominfo['length'][0]*roominfo['height'][0]


    #Now calculate the required HRR that produces flashover
    Q_FO = 378*(1+.021*(surf_area/(roominfo['ventarea'][0]*\
        roominfo['ventheight'][0]**.5)))\
        *roominfo['ventarea'][0]*roominfo['ventheight'][0]**.5

    print Q_FO

    #If/when flashover occurs, the NaN is replaced with a time
    flashover = float('NaN')


    #construct an HRR dictionary holding curves for each item
    hrrdic = {}

    #Incluedes a draw for multiple curves if applicable
    for itemtype in iteminfo.index:
        hrrdic[itemtype] = np.loadtxt(configuse['itemhrrdataloc']\
        +itemtype+ str(np.random.randint(0,iteminfo.numcurves[itemtype])) \
        + '.csv',delimiter=',')

    

    #Retrieving a list of items in the room
    rawlayout= pd.read_csv(configuse['layoutloc'] \
    + '/' + layout + '.csv')
    itemlist = rawlayout.item


    #Then determine which item ignited first
    firstign = False
    for item in range(0,len(itemlist)):
        if (rawlayout['x1'][item] <= fireloc[0] <= rawlayout['x2'][item]) \
        and (rawlayout['y1'][item] <= fireloc[1] <= rawlayout['y2'][item]):
            firstign = item

    #If no item has the fire coordinate, raise an error
    if firstign == False:
        raise fireLocError("The specified fire location is not on an item.")

    num_items = len(itemlist)

    #holds a list of ignition times
    fire_list = np.ones(len(itemlist))*configuse['simtime']
    fire_list[firstign] = 0

    #Holds the location of ignition
    ign_loc = [None]*len(itemlist)
    ign_loc[firstign] = fireloc

    source_list = [None]*len(itemlist)
    source_list[firstign] = fireloc


    #Creating a ist of accumulated FTP for each item
    FTP = np.zeros(num_items) #FTP for each item


    for t, time in enumerate(timelist):

        #The list of incident fluxes. Resets each time step
        incident_flux = np.zeros(num_items)
        #The current HRR of the total fire, resets at each time step
        HRR = 0

        for f,fire in enumerate(itemlist):
            if fire_list[f] != configuse['simtime']:

                #Add all fires' contribution to overall HRR
                #Negative times (unignited) are given zero in interp
                HRR = HRR + np.interp(time - fire_list[f]\
                ,hrrdic[fire][:,0],hrrdic[fire][:,1])


                #Check if flashover conditions are met
                if HRR > Q_FO and np.isnan(flashover):
                    flashover = time
                    #If anything hasn't ignited, it has now
                    fire_list[fire_list==configuse['simtime']] = time
                for i,item in enumerate(itemlist):
                    #Only executes if ith item hasn't ignited 
                    #AND fth item has ignited
                    
                    if fire_list[i] == configuse['simtime']:

                        #Use point source to find incident flux on ith item

                        #Calculate distance from f source to nearest i edge
                        if rawlayout.iloc[i]['x1'] <= source_list[f][0] <= rawlayout.iloc[i]['x2']:
                            x_edge = source_list[f][0]
                        elif rawlayout.iloc[i]['x1'] > source_list[f][0]:
                            x_edge = rawlayout.iloc[i]['x1']
                        elif rawlayout.iloc[i]['x2'] < source_list[f][0]:
                            x_edge = rawlayout.iloc[i]['x2']


                        if rawlayout.iloc[i]['y1'] <= source_list[f][1] <= rawlayout.iloc[i]['y2']:
                            y_edge = source_list[f][1]
                        elif rawlayout.iloc[i]['y1'] > source_list[f][1]:
                             y_edge = rawlayout.iloc[i]['y1']
                        elif rawlayout.iloc[i]['y2'] < source_list[f][1] :
                            y_edge = rawlayout.iloc[i]['y2']
                        dist = np.linalg.norm(np.asarray(source_list[f])-[x_edge,y_edge])


                        incident_flux[i] = incident_flux[i] \
                        + iteminfo.radfrac[fire]\
                        *np.interp(time - fire_list[f]\
                        ,hrrdic[fire][:,0],hrrdic[fire][:,1]) \
                        /(4*np.pi*dist**2)

                        #Then use FTP model
                        if incident_flux[i] > iteminfo.qcrit[item]:
                            FTP[i] = FTP[i] + (incident_flux[i]\
                            -iteminfo.qcrit[item])\
                            **iteminfo.n[item]*configuse['timestep']
                 
                            if FTP[i] > iteminfo.FTP[item]:
                                fire_list[i] = time
                                ign_loc[i]= [x_edge,y_edge]
                                source_list[i] = [x_edge,y_edge]

                if np.isnan(flashover):
                    #Finally update source of f
                    xCentroid = (rawlayout['x1'][f] + rawlayout['x2'][f]) / 2.0
                    yCentroid = (rawlayout['y1'][f] + rawlayout['y2'][f]) / 2.0

                    #The value of the max HRR
                    maxHrr = max(hrrdic[fire][:,1])

                    #The index where it occurs
                    indexmax = np.where(hrrdic[fire][:,1]==maxHrr)[0][0]
                    fraclist = hrrdic[fire][:,1]/maxHrr

                    frac = np.interp(time - fire_list[f],\
                        hrrdic[fire][:,0][0:indexmax+1],fraclist[0:indexmax+1])

                    #Uses simple interpolation to migrate source toward centroid
                    xSource = ign_loc[f][0] + \
                    frac*(xCentroid-ign_loc[f][0])

                    ySource = ign_loc[f][1] + \
                    frac*(yCentroid-ign_loc[f][1])

                    source_list[f] = [xSource,ySource]

    pdb.set_trace()
    return [flashover,fire_list]

    
                                           
