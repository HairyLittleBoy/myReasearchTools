from abaqus import mdb
from math import *
from abaqusConstants import *
from caeModules import *
from abaqus import session

# ---------------------------------------------------------------------
# The following scripts are to plot the inverse pole figure contour for
#     the Euler angles obtained from CP-FEM simulation result
# ---------------------------------------------------------------------
# *Step-1
# Extract the Euler angle from odb file of CP-FEM, this script should be run in the abaqus
jobName = 'E:\MyPapers\TiSpherInden\CPFEMSimu\Job-1'
odb = session.openOdb(name=jobName + '.odb')
EAFile = open('EAofGP.csv', 'w')
EA_phi1 = odb.steps['Step-1'].frames[-1].fieldOutputs['SDV130'].values
EA_phi = odb.steps['Step-1'].frames[-1].fieldOutputs['SDV131'].values
EA_phi2 = odb.steps['Step-1'].frames[-1].fieldOutputs['SDV132'].values
for i in range(len(EA_phi1) + 1):
    if i == 0:
        EAFile.write('%2s, %2s, %5s, %5s, %5s\n'
                     % ('Element', 'GPoint', 'phi1', 'PHI', 'phi2'))
    else:
        EAFile.write('%6i, %6i, %8.4F, %8.4F, %8.4F\n'
                     % (EA_phi1[i - 1].elementLabel, EA_phi1[i - 1].integrationPoint,
                        EA_phi1[i - 1].data, EA_phi[i - 1].data, EA_phi2[i - 1].data))
else:
    EAFile.close()
# Extract the coordinates of integration Pt from odb file of CP-FEM,
# this script should be run in the abaqus
# It should be noted that abaqus does not
# output the coordinate of int Pt by default and there is NO GUI access
# to set the output of coordinate of int Pt. To output the coordinate of int Pt,
# .inp file should be modified:
# <*Element Output, directions=YES
# LE, PE, PEEQ, PEMAG, S, SDV,COORD>
# i.e. add "COORD" under "*Element Output"
COORDintPtFile = open('COORDofGP.csv', 'w')
COORDIntPt = odb.steps['Step-1'].frames[-1].fieldOutputs['COORD'].values
for i in range(len(COORDIntPt) + 1):
    if i == 0:
        COORDintPtFile.write('%2s, %2s, %5s, %5s, %5s\n' % ('Element', 'GPoint', 'coordX', 'coordY', 'coordZ'))
    else:
        COORDintPtFile.write(
            '%6i, %6i, %8.4F, %8.4F, %8.4F\n' %
            (COORDIntPt[i - 1].elementLabel, COORDIntPt[i - 1].integrationPoint,
             COORDIntPt[i - 1].data[0], COORDIntPt[i - 1].data[1],
             COORDIntPt[i - 1].data[2]))
else:
    COORDintPtFile.close()

# *Step-2
# Coordinate on inverse pole figure triangular corresponded by the Euler angle
# This script extract result from odb file and draw the 3d contour of the orientation of grain according to the color
#      on the inverse pole figure
# This script needs to be run within the virtual machine because the DAMASK package cannot be installed in the windous
#      platform
import numpy as np
import damask
import pandas as pd

def stereographicProjection(direction):
    d_ = direction / np.linalg.norm(direction)
    return d_[:2] / (1 + abs(d_[2]))

df = pd.read_csv('E:\MyPapers\TiSpherInden\CPFEMSimu\EAofGP.csv')
crystal_direction = np.array([1, 0, 0], dtype=float)
for i in range(len(df.values)):
    Eulers = [df.values[i, 2], df.values[i, 3], df.values[i, 4]]
    r_lab_to_cryst = damask.Rotation.from_Euler_angles(Eulers, degrees=True)
    PF_coords = stereographicProjection(r_lab_to_cryst.inversed() * crystal_direction)


# -------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
df = pd.read_csv('E:\MyPapers\TiSpherInden\CPFEMSimu\EAofGP.csv')
gp = pd.read_csv('E:\MyPapers\TiSpherInden\CPFEMSimu\COORDofGP.csv')
fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(111,projection='3d')
color = np.array([df.values[:, 2],df.values[:, 3],df.values[:, 4]]).T
for i in range(len(df.values)):
     ax.scatter(gp.values[i, 2],gp.values[i, 3],gp.values[i, 4],s=10,c=color[1,:]/255.0)
plt.show


session.Spectrum(name="custom",colors=('#000004','#010107','#02020C','#030312','#050417','#07051D','#0A0722',
                                           '#0D0828','#100A2E','#130B34','#170B3B','#1B0C41','#1F0C47','#230C4D',
                                           '#270B52','#2C0B57','#300A5C','#350A60','#390963','#3E0965','#420A68',
                                           '#460B69','#4A0C6B','#4E0D6C','#530E6D','#57106D','#5B116E','#5F136E',
                                           '#63146E','#67166E','#6B176E','#6F196E','#731A6E','#771C6D','#7B1D6D',
                                           '#7F1F6C','#83206B','#87216A','#8B2369','#8F2468','#932667','#982766',
                                           '#9C2964','#A02A63','#A42C61','#A82D5F','#AB2F5D','#AF315B','#B33359',
                                           '#B73557','#BB3755','#BF3952','#C23B50','#C63D4D','#C9404B','#CD4248',
                                           '#D04545','#D34842','#D74B3F','#DA4E3D','#DD513A','#DF5437','#E25734',
                                           '#E55B31','#E75F2D','#E9622A','#EC6627','#EE6A24','#F06E21','#F1721D',
                                           '#F3771A','#F57B17','#F67F13','#F78410','#F8880D','#F98D0A','#FA9207',
                                           '#FB9606','#FB9B06','#FCA007','#FCA50A','#FCA90E','#FCAE13','#FCB318',
                                           '#FBB81D','#FBBD23','#FAC229','#F9C72F','#F8CC36','#F7D13D','#F6D644',
                                           '#F5DB4C','#F4E055','#F3E55E','#F2E968','#F1EE72','#F2F27C','#F3F587',
                                           '#F5F991','#F8FC9B','#FCFFA4'))
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(spectrum='custom')



