function matrix = euler2Matrix(eulerAngle)
% eulerAngle in degree
 phi1 = eulerAngle(1);
 phi = eulerAngle(2);
 phi2 = eulerAngle(3);
 matrix(1,1) = cosd(phi1)*cosd(phi2) - sind(phi1)*sind(phi2)*cosd(phi);
 matrix(1,2) = sind(phi1)*cosd(phi2) + cosd(phi1)*sind(phi2)*cosd(phi);
 matrix(1,3) = sind(phi2)*sind(phi);
 matrix(2,1) = -cosd(phi1)*sind(phi2)-sind(phi1)*cosd(phi2)*cosd(phi);
 matrix(2,2) = -sind(phi1)*sind(phi2)+cosd(phi1)*cosd(phi2)*cosd(phi);
 matrix(2,3) = cosd(phi2)*sind(phi);
 matrix(3,1) = sind(phi1)*sind(phi);
 matrix(3,2) = -cosd(phi1)*sind(phi);
 matrix(3,3) = cosd(phi);
 

end

