function [stress_fail,shear_fail] = failure()

%determine shear fail and stress fail based on sample data
%define young's modulus 
E_balsa = 3.2953; %[GPa]
E_foam = .035483; %[Gpa]

%read data
data = xlsread('TestData.xlsx');

%break up data for when member broke in shear or moment
MBreak = zeros(14,5);
VBreak = zeros(5);
c = 1;
%[3 4 8 9 12 13 15 17 18 19 21 22 23 24]
for i = [3 4 8 9 12 13 15 17 18 19 21 22 23 24]
    MBreak(c,:) = data(i,:);
    c = c+1;
end
c = 1;
for i = [1 2 5 14 20]
    VBreak(c,:) = data(i,:);
    c = c+1;
end

%determine the moments for each test
for i = 1:14
    F = MBreak(i,2);
    M(i) = (1/2)*F*MBreak(i,3);
end
a = MBreak(:,3);
w = MBreak(:,4);
%exclude outliers
for i = 1:14
    if M(i) > 12 || M(i) < 7
        M(i) = 0;
    end
end
p = find(M ~= 0);
%determine a, w, and max moment values
M = M(p);
a = a(p);
w = w(p);

I_foam = zeros(1,length(M));
I_balsa = zeros(1,length(M));

%determine moment of inertia for balsa wood and foam
for i = 1:length(M)
    I_foam(i) = (w(i)*0.01905^3)/12;
    I_balsa(i) = (w(i)*0.00211667^3)/12 + w(i)*0.00211667*.009923^2;
end

stress_fail = zeros(1,length(M));
%Use flexure formula to calculate mass stress
for i = 1:length(M)
    stress_fail(i) = -M(i)/(I_balsa(i)+I_foam(i)*(E_foam/E_balsa));
end

V = zeros(1,length(VBreak));
%determine max shear
shear_fail = zeros(1,length(VBreak));
for i = 1:length(VBreak) 
    V(i) = VBreak(i,2)/2;
    shear_fail(i) = (3/2)*(V(i)/(VBreak(i,4)*0.01905)); %area in meters^2, shear in newtons
end
%orginize from least to greatest
shear_fail = sortrows(shear_fail',1);
shear_fail = shear_fail';
stress_fail = sortrows(abs(stress_fail)',1);
stress_fail = stress_fail';
%note: V and M are max shear and moments