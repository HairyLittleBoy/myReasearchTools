function schmidFactor = schmidCalc(slpDir, slpPln, strsDir)
% slpDir, slpPln, strsDir are all in cartesian coordinate system
lamda = dot(slpDir,strsDir)/(norm(slpDir,2)*norm(strsDir,2));
phi = dot(slpPln,strsDir)/(norm(slpPln,2)*norm(strsDir,2));
schmidFactor = lamda*phi;
end