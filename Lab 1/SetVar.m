%ASEN 2001 Statics, Structures, and Materials
%ReadFile: define matrices determined by text file sections
%Ben Hagenau
%Made: 8/29/16

function [tracker,numlines,numF,numM,coorF,magdirF,coorM,magdirM,locS,FSdir,MSdir,text] = SetVar(text) %input text
numFM = str2num(text{1});%set the first cell, representing number of F and M, as an array
%determine number of useful rows in the file
numlines=length(text);
%create row number tracker
tracker=2;
%determine number of forces and moments individually
numF=numFM(1);
numM=numFM(2);
%Determine the number of lines each section will have and turn them into their own matrices
%Coordinates of the points at which external forces are applied: number of forces = number of lines
for i=tracker:tracker+numF-1
    coorF(i,:)=str2num(text{i});
    tracker=tracker+1;
end 
%Magnitude and direction of external forces: number of lines = number of forces
for i=tracker:tracker+numF-1
       magdirF(i,:)=str2num(text{i});
       tracker=tracker+1;
end
%Location of external couple moments: number of lines = number of moments
for i=tracker:tracker+numM-1
    coorM(i,:)=str2num(text{i});
    tracker=tracker+1;
end
%magnitude and direction of couple moments: number of lines = number of moments
for i=tracker:tracker+numM-1
    magdirM(i,:)=str2num(text{i});
    tracker=tracker+1;
end
%Number of supports: remaining number of lines /2
lines_left=(19-tracker);
for i=tracker:tracker+lines_left/2
    locS(i,:)=str2num(text{i});
    tracker=tracker+1;
end

fcount=0;
mcount=0;
%Type and direction of support: number of lines = to number of supports (edit letters out of string)
%Determine which are F and which are M
for i=tracker:numlines
    if text{i}(1)=='F'
        fcount=fcount+1;
    else
        mcount=mcount+1;
    end
end
%get rid of all characters in the last set lines
for i=tracker:tracker+lines_left/2
    text{i}(1)=[];
end
%create 2 vectors for F and M supports
for i=tracker:tracker+fcount-1
    FSdir(i,:)=str2num(text{i});
end
for i=numlines-mcount+1:numlines
    MSdir(i,:)=str2num(text{i});
end

%Clear full zero (unecessary) rows from all matrices 

    