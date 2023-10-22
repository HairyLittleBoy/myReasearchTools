function vectorcartsn = hex2cartsnDirYP(vectorhex)
% 电子背散射衍射技术及其应用(杨平) 4.1
caRatio = 1.587;
u = vectorhex(1);
v = vectorhex(2);
t = vectorhex(3);
w = vectorhex(4);
vectorcartsn(1) = sqrt(3)/2*(2*u+v);
vectorcartsn(2) = 3/2*v;
vectorcartsn(3) = w*caRatio;
vectorcartsn = vectorcartsn;
%vectorcartsn = vectorcartsn./norm(vectorcartsn,2);
end

