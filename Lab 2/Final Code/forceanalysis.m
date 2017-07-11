function [barforces,reacforces,jointweight]=forceanalysis(joints,connectivity,reacjoints,reacvecs)
% function [barforces,reacforces]=forceanalysis(joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
%
% compute forces in bars and reaction forces
%
% input:  joints       - coordinates of joints
%         connectivity - connectivity 
%         reacjoints   - joint id where reaction acts on
%         reacvecs     - unit vector associated with reaction force
%         loadjoints   - joint id where external load acts on
%         loadvecs     - load vector
%
% output: barforces    - force magnitude in bars
%         reacforces   - reaction forces
%
% Author: Kurt Maute, Sept 21 2011

% extract number of joints, bars, reactions, and loads
numjoints = size(joints,1);
numbars   = size(connectivity,1);
numreact  = size(reacjoints,1);

% number of equations
numeqns = 3 * numjoints;

% allocate arrays for linear system
Amat = zeros(numeqns);
bvec = zeros(numeqns,1);

% build Amat - loop over all joints

for i=1:numjoints
    
   % equation id numbers
   idx = 3*i-2;
   idy = 3*i-1;
   idz = 3*i;
   
   % get all bars connected to joint
   [ibar,ijt]=find(connectivity==i);
   
   % loop over all bars connected to joint
   for ib=1:length(ibar)
       
       % get bar id
       barid=ibar(ib);
       
       % get coordinates for joints "i" and "j" of bar "barid"
       joint_i = joints(i,:);
       if ijt(ib) == 1
           jid = connectivity(barid,2);
       else
           jid = connectivity(barid,1);
       end
       joint_j = joints(jid,:);
       
       % compute unit vector pointing away from joint i
       vec_ij = joint_j - joint_i;
       uvec   = vec_ij/norm(vec_ij);
       
       % add unit vector into Amat
       Amat([idx idy idz],barid)=uvec;
   end
end

% build contribution of support reactions 
for i=1:numreact
    
    % get joint id at which reaction force acts
    jid=reacjoints(i);

    % equation id numbers
    idx = 3*jid-2;
    idy = 3*jid-1;
    idz = 3*jid;

    % add unit vector into Amat
    Amat([idx idy idz],numbars+i)=reacvecs(i,:);
end

%find bar info
barinfo=zeros(numbars,3);
for i=1:numbars
    barinfo(i,1)=i;
    
    P1=connectivity(i,1);
    P2=connectivity(i,2);
    
    Coord1=joints(P1,1:3);
    Coord2=joints(P2,1:3);
    
    Length=((Coord2(1)-Coord1(1)).^2+((Coord2(2)-Coord1(2)).^2 +((Coord2(3)-Coord1(3)).^2))).^(1/2);
    barinfo(i,2)=Length;
    Length=Length-2*.01272; %length = length - 2*diameterball 
    
    barinfo(i,3) = 0.0339*Length; %use kg/m density
    barinfo(i,3) = barinfo(i,3) + 2*0.001679;
end

% build load vector
jointweight =zeros(numjoints,1);
for i=1:numjoints
    for j=1:numbars
        if connectivity(j,1)==i
            jointweight(i)=jointweight(i)+(barinfo(j,3)/2);
        elseif connectivity(j,2)==i
            jointweight(i)=jointweight(i)+(barinfo(j,3)/2);
        end
    end
    jointweight(i)=jointweight(i)+0.008358; %weight of magnet in kg
end

jointweight=jointweight*9.807;
for i=1:numjoints
    bvec(3*i)=jointweight(i);
end

% check for invertability of Amat
if rank(Amat) ~= numeqns
    error('Amat is rank defficient: %d < %d\n',rank(Amat),numeqns);
end

% solve system
xvec=Amat\bvec;

% extract forces in bars and reaction forces
barforces=xvec(1:numbars);
reacforces=xvec(numbars+1:end);

end