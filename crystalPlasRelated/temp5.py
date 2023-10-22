for i in range(125):
    mdb.models['n125-id1'].Material(name='Material-'+str(i+1))
    mdb.models['n125-id1'].materials['Material-'+str(i+1)].UserMaterial(mechanicalConstants=(0.0, ))
    mdb.models['n125-id1'].materials['Material-'+str(i+1)].Depvar(n=200)


    mdb.models['n125-id1'].HomogeneousSolidSection(name='Section-'+str(i+1),material='Material-'+str(i+1), thickness=None)
    region = mdb.models['n125-id1'].parts['TESS'].sets['POLY'+str(i+1)]
    p = mdb.models['n125-id1'].parts['TESS']
    p.SectionAssignment(region=region, sectionName='Section-'+str(i+1))