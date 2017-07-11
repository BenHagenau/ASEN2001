function truss3dweight(inputfile,outputfile)
% function truss2d(inputfile,outputfile)
%
% Analysis of 2-D statically determinate truss
%
% Input:  inputfile  - name of input file
%         outputfile - name of output file
%
% Author: Kurt Maute for ASEN 2001, Sept 21 2011

% read input file
[joints,connectivity,reacjoints,reacvecs]=readinput(inputfile);

% compute forces in bars and reactions
[barforces,reacforces,loadjoints]=forceanalysis(joints,connectivity,reacjoints,reacvecs);

% write outputfile
writeoutput(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints);

% plot truss (used in Lab 2)
joints3D=zeros(size(joints,1),3);
joints3D(:,1:3)=joints;
plottruss(joints3D,connectivity,barforces,reacjoints,(1/3)*[0.025,0.04,0.05],[1 1 0 0])

end
