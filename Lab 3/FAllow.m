%fucntion to determine the allowable stress or shear

function [F_allow] = FAllow(FS,F_fail) 
F_allow = (F_fail/FS);