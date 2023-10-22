function vectorcartsn = hex2cartsnDir(vectorhex)
caRatio = 1.587;
r = sqrt(vectorhex(1)^2+vectorhex(2)^2-vectorhex(1)*vectorhex(2));
theta = atand(sqrt(3)*vectorhex(2)/(2*vectorhex(1)-vectorhex(2)));
if (vectorhex(1)>0&&vectorhex(2)>=0)
    theta = theta;
elseif (vectorhex(1)<=0&&vectorhex(2)>=0)
    theta = 180+theta;
elseif (vectorhex(1)<0&&vectorhex(2)<0)
    if abs(vectorhex(1))>=abs(vectorhex(2))
        theta = 180+theta;
    elseif abs(vectorhex(1))<abs(vectorhex(2))
        theta = 360+theta;
    end
elseif (vectorhex(1)>=0&&vectorhex(2)<0)
    theta = 360+theta;
end
vectorcartsn(1) = r*cosd(theta);
vectorcartsn(2) = r*sind(theta);
vectorcartsn(3) = vectorhex(3)*caRatio;
vectorcartsn = vectorcartsn./norm(vectorcartsn,2);
end

