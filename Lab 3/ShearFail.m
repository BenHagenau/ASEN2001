%determine the failure shear
function [shear_fail] = ShearFail(F,w)
%calculate failing shear
shear_fail = (3/2)*((F/2)/(w*0.01905));