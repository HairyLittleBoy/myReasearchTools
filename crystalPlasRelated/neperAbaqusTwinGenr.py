# the twinned grain number
# (GOOD GRAINS) POLY3,POLY4,POLY5,POLY9,POLY10,POLY20,POLY29,POLY52,POLY61,POLY85,POLY86,POLY99,POLY114,POLY125


from abaqus import *
from abaqusConstants import *
import numpy as np

def neperAbaqusTwinGenr(grainNum, twinLamellaThickness, twinPlaneNormal,
                        twinlength, twinStartPoint, twinSetName) :
    # twinLamellaThickness : mm
    # twinLamellaThickness : t
    # twinPlaneNormalUnit : [nx, ny, nz]
    # twinStartPoint : [xs, ys, zs]
    # twinStartPoint_1 : [x1, y1, z1]
    # twinStartPoint_2 : [x2, y2, z2]
    # how to select the elements of twin :
    # the elements between plane_1 and plane_2 will be selected,
    #      plane_1's normal is [nx, ny, nz] and passes twinStartPoint_1 [x1, y1, z1]
    #      plane_2's normal is [nx, ny, nz] and passes twinStartPoint_2 [x2, y2, z2]
    #      plane_1's equation is nx*X + ny*Y + nz*Z + D1 = 0, where D1 = -nx*x1 - ny*y1 -nz*z1
    #      plane_2's equation is nx*X + ny*Y + nz*Z + D2 = 0, where D2 = -nx*x2 - ny*y2 -nz*z2

    twinPlaneNormalUnit = twinPlaneNormal / np.linalg.norm(twinPlaneNormal)
    twinStartPoint_1 = 0.5 * twinPlaneNormalUnit * twinLamellaThickness + twinStartPoint
    twinStartPoint_2 = -0.5 * twinPlaneNormalUnit * twinLamellaThickness + twinStartPoint

#    D1 = twinPlaneNormalUnit * twinStartPoint_1
#    D2 = twinPlaneNormalUnit * twinStartPoint_2

    twinedGrain = mdb.models['n125-id1'].parts['TESS'].sets[grainNum]
    eleNumIntwinedGrain = len(twinedGrain.elements)

    eleInTwin = []
    distanceNodeToPlane_1 = np.zeros((eleNumIntwinedGrain,8))
    distanceNodeToPlane_2 = np.zeros((eleNumIntwinedGrain,8))
    for i in range(eleNumIntwinedGrain):
        for j in range(8):
            nodeCoord = twinedGrain.elements[i].getNodes()[j].coordinates
            distanceNodeToPlane_1[i,j] = abs(np.dot([nodeCoord - twinStartPoint_1], twinPlaneNormalUnit)[0])
            distanceNodeToPlane_2[i,j] = abs(np.dot([nodeCoord - twinStartPoint_2], twinPlaneNormalUnit)[0])
            distanceNodeToStartPoint = np.linalg.norm([nodeCoord - np.array(twinStartPoint)])
            distanceNodeToPlanes = distanceNodeToPlane_1[i,j] + distanceNodeToPlane_2[i,j]
            if abs(distanceNodeToPlanes-twinLamellaThickness) < 1.e-4 and distanceNodeToStartPoint < twinlength :
                eleInTwin.append(twinedGrain.elements[i].label)
    eleInTwin = np.unique(eleInTwin)            
    twinSequence = mdb.models['n125-id1'].parts['TESS'].elements.sequenceFromLabels(eleInTwin)
    mdb.models['n125-id1'].parts['TESS'].Set(elements=twinSequence, name=twinSetName)
    return eleInTwin
    
# -----------------------------------------------------------------------------------------------------
twinInfor = {}
twinInfor = {       'grainNum':['POLY5','POLY9','POLY10','POLY20','POLY29','POLY52',
                                'POLY61','POLY85','POLY86','POLY99','POLY114','POLY125'],
              'twinStartPoint':[[740.E-03,1.,400.E-03],
                               [265.06E-03,927.711E-03,674.699E-03],
                               [590.361E-03,1.,0.],
                               [530.12E-03,385.542E-03,349.398E-03],
                               [361.446E-03,783.133E-03,397.59E-03],
                               [409.639E-03,554.217E-03,265.06E-03],
                               [819.277E-03,349.398E-03,421.687E-03],
                               [614.458E-03,734.94E-03,722.892E-03],
                               [771.084E-03,638.554E-03,746.988E-03],
                               [987.952E-03,987.952E-03,216.867E-03],
                               [1.,626.506E-03,144.578E-03],
                               [60.241E-03,783.133E-03,361.446E-03]],
             'twinPlaneNormal':[[-1.,  1.,  0.5],
                               [ 1.,  1., -1.],
                               [-0.5, 1.,  1.],
                               [ 0.1, 1.,  0.1],
                               [ 0.1, 1., -0.1],
                               [ 0.6, 1.5, 1.],
                               [-0.5, 1.,  0.5],
                               [-0.2, 1., -1.],
                               [-1.,  1., -1.],
                               [-1.,  1.,  0.5],
                               [-1.,  1.,  1.],
                               [ 1.,  1.,  1.]]}


twinLamellaThickness = 0.03
twinlength = 0.5
totalEleNumInTwin = []

for i in range(len(twinInfor['grainNum'])):

    eleInTwin = []
    eleInTwin = neperAbaqusTwinGenr(twinInfor['grainNum'][i], twinLamellaThickness, 
                                    twinInfor['twinPlaneNormal'][i], twinlength, 
                                    twinInfor['twinStartPoint'][i], 'twin'+twinInfor['grainNum'][i])

    eleInTwinFile = open('eleInTwin'+twinInfor['grainNum'][i]+'.txt', 'w')
    for j in range(len(eleInTwin) + 1):
        if j == 0:
            eleInTwinFile.write('%8s,%5i\n' % ('(/', eleInTwin[j]))
        elif j % 9 == 0 and j != 0 and j != len(eleInTwin):
            eleInTwinFile.write('\n %5s, %5i' % ('&', eleInTwin[j]))
        elif j == len(eleInTwin):
            eleInTwinFile.write('%2s' % ('/)'))
        else:
            eleInTwinFile.write('%1s%5i' % (',', eleInTwin[j]))
    totalEleNumInTwin.append(len(eleInTwin))
    
    eleInTwinFile.close()
    
    
    