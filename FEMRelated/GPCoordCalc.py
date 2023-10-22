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
#   nodesCoord : [[x0,y0],[x1,y1] ... [x7,y7]]
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



import numpy as np
import math
elemLabel = 1
nodesCoord = [np.array([0., 0., 0.0]),
              np.array([10., 0., 0.0]),
              np.array([10., 10., 0.0]),
              np.array([0., 10., 0.0]),
              np.array([5.0, 0., 0.0]),
              np.array([10.0, 5., 0.0]),
              np.array([5.0, 10., 0.0]),
              np.array([0., 5., 0.0])]

nodesCoord = [np.array([0.000341497594490647, -8.40511226654053, 0.0]),
              np.array([0.000695057562552392, -8.40511226654053, 0.0]),
              np.array([0.000695057562552392, -8.11802768707275, 0.0]),
              np.array([0.000341497594490647, -8.11802768707275, 0.0]),
              np.array([0.0, -8.55372524261475, 0.0]),
              np.array([0.000170748797245324, -8.11802768707275, 0.0]),
              np.array([0.0, -8.26156997680664, 0.0]),
              np.array([0.000695057562552392, -8.26156997680664, 0.0])]

GPPoints = [-1./math.sqrt(3.), 1./math.sqrt(3.)]
weights = [1.,1.]
GPCoordsElem = GPCoodCalc(1,nodesCoordOfElem)