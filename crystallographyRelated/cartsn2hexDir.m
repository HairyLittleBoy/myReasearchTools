function vextorHex4 = cartsn2hexDir(vextorCartsn)
caRatio = 1.587;
vextorHex4(4) = vextorCartsn(3)/ caRatio;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
x0 = [0.5,0.5];
lb = [-1,-1];
ub = [1,1];
[x,~,~,~,~,~,~] = ...
fmincon(@(x)ErrorCartsn2hexDir(x,vextorCartsn),x0,[],[],[],[],lb,ub,[],options);
vextorHex3(1) = x(1);
vextorHex3(2) = x(1);
HCPIdx ([vextorHex3(1),vextorHex3(2),vextorHex4(4)], 'dir')

end

