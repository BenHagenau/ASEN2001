clear all
clc
close all

%Youngs Modulus
E_balsa = 3.2953; %[GPa]
E_foam = .035483; %[Gpa]
%define factor of safety
FS = 1.7;
%determine shear fail and stress fail
[stress_fail,shear_fail] = failure();
%define allowable shear of stress (use greatest value of shear fail and
%moment fail)
[stress_allow] = FAllow(FS,abs(stress_fail(end)));
[shear_allow] = FAllow(FS,abs(shear_fail(end)));

%determine maximum pressure throughout the pressure field (applied at
%middle)
%4 inchs = 0.1016 m
%p0 = stress_allow*0.00235974;

%Determine p0
p0 = (stress_allow/(FS*.0667))*10^-6;
%p0 = 4.8455;
%I as a function of width
I_f = @(w) (w*0.01905^3)/12;
I_b = @(w) (w*0.00211667^3)/12 + w*0.00211667*.009923^2;

%determine width values
c = 1;
x = linspace(-0.4572,0.4572,100);
for i = x
    wV(c) = -(((.0106*p0*sqrt(1 - 4.784*i^2) - .00354*p0*(1 - 4.784*i.^2).^(3/2)...
       + .0232*p0*i*asin(2.187*i))*.9144)/(stress_allow*(((1/16).^3)/12)...
       + (E_balsa/E_foam)*((.75^3)/12)+0.00211667*.009923^2));
    c = c+1;
    %wM = 
end


%plot
figure
plot(x,wV,'o')
hold on
hold off