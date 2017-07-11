%calculate the failure stress

function [stress_fail] = StressFail(F,w,a)
%define young's modulus 
E_balsa = 3.2953; %[GPa]
E_foam = .035483; %[Gpa]

%calculate moments of inertia
I_foam = (w*0.01905^3)/12;
I_balsa = (w*0.00211667^3)/12 + w*0.00211667*.009923^2;

%calculate moment
M = (1/2)*F*a;

%calculate failing stress
stress_fail = -M/(I_balsa+I_foam*(E_foam/E_balsa));

