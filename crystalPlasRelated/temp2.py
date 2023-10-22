import scipy.io as sio
import numpy as np

matfn='E:/MyPapers/TiSpherInden/CPFEMSimu/EulerAngles125Grains.mat' #read.mat
data = sio.loadmat(matfn)

# load_h=data['h']                     #read the h array in .mat,the h in Miller index
# load_k=data['k']
# load_l=data['l']

# load_x=data['x']
# load_y=data['y']
# load_z=data['z']


for i in range(125):
# the pattern in abaqus for defining the material parameters
# EulerAngle = (/props(158),props(159),props(160)/)
    mechanicalConstants=(162400.0,92000.0,162400.0,66000.0,66000.0,
     180900.0,46700.0,35200.0,46700.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,0.0, 1.0, 0.0,
     0.0, 1.0, 0.0, 0.0, 0.0, 50.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 50.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50.0, 0.001, 0.0,
     0.0, 0.0, 0.0, 0.0,0.0, 200.0, 568.0, 127.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0, 1500.0, 96.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     300.0, 400.0, 240.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
     50.0, 1.0e-05, 0.0, 0.0, 0.0, 0.0, 0.0)

    Constants = list(mechanicalConstants)

    load_euler1 = data['EulerAngles'][i][0]
    load_euler2 = data['EulerAngles'][i][1]
    load_euler3 = data['EulerAngles'][i][2]
#modify the grain orientation at corresponding position    
    Constants[157] = load_euler1
    Constants[158] = load_euler2
    Constants[159] = load_euler3

    mc = tuple(Constants)

    f = open('pureTi125G.txt','a')

    f.write('\n')

    k = str(i+1)

    f.write('mdb.models['+'\'n125-id1\''+'].materials['+'\'Material-'+k+'\']'+'.userMaterial.setValues(mechanicalConstants=')

    for j in range(160):

        print(j)

        if j==0:

           f.write('(')

        var=str(Constants[j])

        f.writelines([var,','])

        if j==159:

            f.write('))')

    f.close()