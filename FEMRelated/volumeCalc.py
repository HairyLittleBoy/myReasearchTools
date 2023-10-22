def GPCoodCalc(elemLabel,nodesCoord):
    import math
#   the (g,h) coordinates of gauss integral points
#   (1) : -sqrt(3)/3,-sqrt(3)/3
#   (2) : sqrt(3)/3,-sqrt(3)/3
#   (3) : sqrt(3)/3,sqrt(3)/3
#   (4) : -sqrt(3)/3,sqrt(3)/3
#                       ^h
#                       |
#            -----------|-----------
#           |           |           |
#           |      *(4) |    *(3)   |
#           |           |           |
#   --------|-----------|-----------|----------->g
#           |           |           |
#           |      *(1) |    *(2)   |
#           |           |           |
#            -----------|-----------
#                       |
#                       |    
   
    GPgh = [[-math.sqrt(3)/3,-math.sqrt(3)/3],[math.sqrt(3)/3,-math.sqrt(3)/3],[math.sqrt(3)/3,math.sqrt(3)/3],[-math.sqrt(3)/3,math.sqrt(3)/3]]
    GPCoord = []
    for i in range(4):
        g = GPgh[i][0]
        h = GPgh[i][1]

        GPCoordx = -0.25*(1-g)*(1-h)*(1+g+h)*nodesCoord[0][0] - 0.25*(1+g)*(1-h)*(1-g+h)*nodesCoord[1][0] \
                       -0.25*(1+g)*(1+h)*(1-g-h)*nodesCoord[2][0] - 0.25*(1-g)*(1+h)*(1+g-h)*nodesCoord[3][0] \
                       +0.5*(1-g)*(1+g)*(1-h)*nodesCoord[4][0] + 0.5*(1-h)*(1+h)*(1+g)*nodesCoord[5][0] \
                       +0.5*(1-g)*(1+g)*(1+h)*nodesCoord[6][0] + 0.5*(1-h)*(1+h)*(1-g)*nodesCoord[7][0]

        GPCoordy = -0.25*(1-g)*(1-h)*(1+g+h)*nodesCoord[0][1] - 0.25*(1+g)*(1-h)*(1-g+h)*nodesCoord[1][1]\
                       -0.25*(1+g)*(1+h)*(1-g-h)*nodesCoord[2][1] - 0.25*(1-g)*(1+h)*(1+g-h)*nodesCoord[3][1]\
                       +0.5*(1-g)*(1+g)*(1-h)*nodesCoord[4][1] + 0.5*(1-h)*(1+h)*(1+g)*nodesCoord[5][1]\
                       +0.5*(1-g)*(1+g)*(1+h)*nodesCoord[6][1] + 0.5*(1-h)*(1+h)*(1-g)*nodesCoord[7][1]
        GPCoord.append([GPCoordx,GPCoordy])
    GPCoordsElem = {"elementLabel":elemLabel,"GPCoordinates":GPCoord}

    return GPCoordsElem

def elemArea(elemLabel,nodesCoord, GPPoints, weights):
    import math
#   the (g,h) coordinates of gauss integral points
#   (1) : -sqrt(3)/3,-sqrt(3)/3
#   (2) : sqrt(3)/3,-sqrt(3)/3
#   (3) : sqrt(3)/3,sqrt(3)/3
#   (4) : -sqrt(3)/3,sqrt(3)/3
#                       ^h
#                       |
#            -----------|-----------
#           |           |           |
#           |      *(4) |    *(3)   |
#           |           |           |
#   --------|-----------|-----------|----------->g
#           |           |           |
#           |      *(1) |    *(2)   |
#           |           |           |
#            -----------|-----------
#                       |
#                       |
    #   nodesCoord : [[x0,y0],[x1,y1] ... [x7,y7]]
    #   dNdg : [dN1dg, dN2dg ... dN8dg]
    # reference :http://ovainola.github.io/jekyll/update/2016/07/22/Calculating-element-surface-area.html
    area = 0.
    for i in range(len(GPPoints)):
        for j in range(len(GPPoints)):
            g = GPPoints[i]
            h = GPPoints[j]
            wti = weights[i]
            wtj = weights[j]
            dNdg = [-0.25*((-1.)*(1.-h)*(1.+g+h)+(1.-g)*(1.-h)),-0.25*((1.-h)*(1.-g+h)+(1.+g)*(1.-h)*(-1.)),
                    -0.25*((1.+h)*(1.-g-h)+(1.+g)*(1.+h)*(-1.)),-0.25*((-1.)*(1.+h)*(1.+g-h)+(1.-g)*(1.+h)),
                    0.5*((-1.)*(1.+g)*(1.-h)+(1.-g)*(1.-h)),0.5*((1.-h)*(1.+h)),
                    0.5*((-1.)*(1.+g)*(1.+h)+(1.-g)*(1.+h)),0.5*((-1.)*(1.-h)*(1.+h))]
            dNdh = [-0.25*((1.-g)*(-1.)*(1.+g+h)+(1.-g)*(1.-h)),-0.25*((1.+g)*(-1.)*(1.-g+h)+(1.+g)*(1.-h)),
                    -0.25*((1.+g)*(1.-g-h)+(1.+g)*(1.+h)*(-1.)),-0.25*((1.-g)*(1.+g-h)+(1.-g)*(1.+h)*(-1.)),
                    0.5*((1.-g)*(1.+g)*(-1.)),0.5*((-1.)*(1.+h)*(1.+g)+(1.-h)*(1.+g)),
                    0.5*((1.-g)*(1.+g)),0.5*((-1.)*(1.+h)*(1.-g)+(1.-h)*(1.-g))]
            dxdg = 0.
            dxdh = 0.
            dydg = 0.
            dydh = 0.
            for k in range(8):
                dxdg = dxdg + dNdg[k]*nodesCoord[k][0]
                dydg = dydg + dNdg[k]*nodesCoord[k][1]
                dxdh = dxdh + dNdh[k]*nodesCoord[k][0]
                dydh = dydh + dNdh[k]*nodesCoord[k][1]
            detJacobi = dxdg * dydh - dydh * dxdh
            area = area + wti * wtj * detJacobi
    AreaElem = {"elementLabel": elemLabel, "Area": area}
    return AreaElem


def elemVolume(elemLabel,nodesCoord, GPPoints, weights):
    import math
#   the (g,h) coordinates of gauss integral points
#   (1) : -sqrt(3)/3,-sqrt(3)/3
#   (2) : sqrt(3)/3,-sqrt(3)/3
#   (3) : sqrt(3)/3,sqrt(3)/3
#   (4) : -sqrt(3)/3,sqrt(3)/3
#                       ^h
#                       |
#            -----------|-----------
#           |           |           |
#           |      *(4) |    *(3)   |
#           |           |           |
#   --------|-----------|-----------|----------->g
#           |           |           |
#           |      *(1) |    *(2)   |
#           |           |           |
#            -----------|-----------
#                       |
#                       |
    #   nodesCoord : [[x0,y0],[x1,y1] ... [x7,y7]]
    #   dNdg : [dN1dg, dN2dg ... dN8dg]

    area = 0.
    volume = 0.
    for i in range(len(GPPoints)):
        for j in range(len(GPPoints)):
            g = GPPoints[i]
            h = GPPoints[j]
            wti = weights[i]
            wtj = weights[j]
            dNdg = [-0.25*((-1.)*(1.-h)*(1.+g+h)+(1.-g)*(1.-h)),-0.25*((1.-h)*(1.-g+h)+(1.+g)*(1.-h)*(-1.)),
                    -0.25*((1.+h)*(1.-g-h)+(1.+g)*(1.+h)*(-1.)),-0.25*((-1.)*(1.+h)*(1.+g-h)+(1.-g)*(1.+h)),
                    0.5*((-1.)*(1.+g)*(1.-h)+(1.-g)*(1.-h)),0.5*((1.-h)*(1.+h)),
                    0.5*((-1.)*(1.+g)*(1.+h)+(1.-g)*(1.+h)),0.5*((-1.)*(1.-h)*(1.+h))]
            dNdh = [-0.25*((1.-g)*(-1.)*(1.+g+h)+(1.-g)*(1.-h)),-0.25*((1.+g)*(-1.)*(1.-g+h)+(1.+g)*(1.-h)),
                    -0.25*((1.+g)*(1.-g-h)+(1.+g)*(1.+h)*(-1.)),-0.25*((1.-g)*(1.+g-h)+(1.-g)*(1.+h)*(-1.)),
                    0.5*((1.-g)*(1.+g)*(-1.)),0.5*((-1.)*(1.+h)*(1.+g)+(1.-h)*(1.+g)),
                    0.5*((1.-g)*(1.+g)),0.5*((-1.)*(1.+h)*(1.-g)+(1.-h)*(1.-g))]
            dxdg = 0.
            dxdh = 0.
            dydg = 0.
            dydh = 0.
            for k in range(8):
                dxdg = dxdg + dNdg[k]*nodesCoord[k][0]
                dydg = dydg + dNdg[k]*nodesCoord[k][1]
                dxdh = dxdh + dNdh[k]*nodesCoord[k][0]
                dydh = dydh + dNdh[k]*nodesCoord[k][1]
            detJacobi = dxdg * dydh - dydh * dxdh
            Nixi = -0.25*(1.-g)*(1.-h)*(1.+g+h)*nodesCoord[0][0] - 0.25*(1.+g)*(1.-h)*(1.-g+h)*nodesCoord[1][0]\
                   -0.25*(1.+g)*(1.+h)*(1.-g-h)*nodesCoord[2][0] - 0.25*(1.-g)*(1.+h)*(1.+g-h)*nodesCoord[3][0]\
                   +0.5*(1.-g)*(1.+g)*(1.-h)*nodesCoord[4][0] + 0.5*(1.-h)*(1.+h)*(1.+g)*nodesCoord[5][0]\
                   +0.5*(1.-g)*(1.+g)*(1.+h)*nodesCoord[6][0] +0.5*(1.-h)*(1.+h)*(1.-g)*nodesCoord[7][0]
            volume = volume + 2. * 3.14159265359 * Nixi * detJacobi
    VolumeElem = {"elementLabel": elemLabel, "volume": volume}
    return VolumeElem
    
# ---------------volume of plastic and elastic region in spherical deformation-----------

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import numpy as np
import math

yldStrs = 200.
path = 'H:/abaqusWorkingFiles/WSY'
odbFileName = 'Job-1'
FrameNum = -1

GPPoints = [-1./math.sqrt(3.), 1./math.sqrt(3.)]
weights = [1.,1.]
o1 = session.openOdb(name=path+'/'+odbFileName+'.odb')

nodesAll=o1.rootAssembly.instances['MATERIAL-1'].nodes
elemsAll=o1.rootAssembly.instances['MATERIAL-1'].elements
nodeNumAll = len(nodesAll)
elemNumAll = len(elemsAll)
MisesOfGPAll=session.odbs[o1.name].steps['Step-1'].frames[FrameNum].fieldOutputs['S'].values

MisesOfGPMat = [[0.]*5 for i in range(elemNumAll)]
for i in range(len(MisesOfGPAll)):
    if MisesOfGPAll[i].instance.name == 'MATERIAL-1' :
       MisesOfGPMat[MisesOfGPAll[i].elementLabel-1][MisesOfGPAll[i].integrationPoint-1] = MisesOfGPAll[i].mises
       MisesOfGPMat[MisesOfGPAll[i].elementLabel-1][4] = MisesOfGPAll[i].elementLabel

areaPlastic = 0.
volumePlastic = 0.
areaElastic = 0.
volumeElastic = 0.
elemLabelPlas = []
elemLabelElas = []
for i in range(len(MisesOfGPMat)):
    currentElemLabel = MisesOfGPMat[i][4]
#        print(currentElemLabel)
    nodesCoordOfElem = []
    MisesOfThisElemGP = MisesOfGPMat[currentElemLabel-1][0:4]
    if (np.array(MisesOfThisElemGP)  >= yldStrs - yldStrs * 1.e-2).all() :
        nodesLabelOfThisElem = elemsAll[i].connectivity       
        for j in range(len(nodesLabelOfThisElem)):
            nodesCoordOfElem.append(nodesAll[nodesLabelOfThisElem[j]-1].coordinates)
        areaThisElem = elemArea(currentElemLabel,nodesCoordOfElem, GPPoints, weights)
        volumeThisElem = elemVolume(currentElemLabel,nodesCoordOfElem, GPPoints, weights)
        areaPlastic = areaPlastic + areaThisElem['Area']
        volumePlastic = volumePlastic + volumeThisElem['volume']
#            print(areaThisElem['Area'])
#            print(areaPlastic)
        elemLabelPlas.append(areaThisElem['elementLabel'])
        
    elif (np.array(MisesOfThisElemGP)  < yldStrs - yldStrs * 1.e-1).all():
        nodesLabelOfThisElem = elemsAll[i].connectivity
        for j in range(len(nodesLabelOfThisElem)):
            nodesCoordOfElem.append(nodesAll[nodesLabelOfThisElem[j]-1].coordinates)
        areaThisElem = elemArea(currentElemLabel,nodesCoordOfElem, GPPoints, weights)
        volumeThisElem = elemVolume(currentElemLabel,nodesCoordOfElem, GPPoints, weights)
        areaElastic = areaElastic + areaThisElem['Area'] 
        volumeElastic = volumeElastic + volumeThisElem['volume']        
        elemLabelElas.append(areaThisElem['elementLabel'])   