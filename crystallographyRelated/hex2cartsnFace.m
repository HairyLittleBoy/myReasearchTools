function UnitNormalCartsn = hex2cartsnFace(facehcp4indx)
% this function calculate the unit normal vector in grain cartesian coordinate 
% system of a crystal face expressed by hcp 4 index
caRatio = 1.587;
if facehcp4indx(1)==0
    vector1 = [0,1/facehcp4indx(2),0];
    vector2 = [-1/facehcp4indx(3),-1/facehcp4indx(3),0];
elseif facehcp4indx(2)==0
    vector1 = [1/facehcp4indx(1),0,0];
    vector2 = [-1/facehcp4indx(3),-1/facehcp4indx(3),0];
else
    vector1 = [1/facehcp4indx(1),0,0];
    vector2 = [0,1/facehcp4indx(2),0];
end
if facehcp4indx(4) == 0
    V1_3 = vector1 - vector2;

    V1 = hex2cartsnDir(V1_3);
    V2 = [0,0,1];
else
    vector3 = [0,0,1/facehcp4indx(4)];
    V1_3 = vector1 - vector3;
    V2_3 = vector2 - vector3;
    V1 = hex2cartsnDir(V1_3);
    V2 = hex2cartsnDir(V2_3);
end

NormalCartsn = cross(V1,V2);
if facehcp4indx(1)==0 && facehcp4indx(2)==0
    NormalCartsn = [0,0,caRatio];
end
UnitNormalCartsn = NormalCartsn./norm(NormalCartsn,2);

end

