function UnitNormalCartsn = hex2cartsnFaceYP(facehcp4indx)
% this function calculate the unit normal vector in grain cartesian coordinate 
% system of a crystal face expressed by hcp 4 index
caRatio = 1.587;
h = facehcp4indx(1);
k = facehcp4indx(2);
l = facehcp4indx(4);
NormalCartsn(1) = (2*h+k)*sqrt(3)/3;
NormalCartsn(2) = k;
NormalCartsn(3) = l*(1/caRatio);
UnitNormalCartsn = NormalCartsn;
%UnitNormalCartsn = NormalCartsn./norm(NormalCartsn,2);

end

