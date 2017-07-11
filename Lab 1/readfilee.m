%ASEN 2001 Statics, Structures, and Materials
%ReadFile: inport file and read file
%Ben Hagenau
%Made: 8/29/16

%Define function Fnumber,Mnumber,Fmagdir,Mmagdir,Slocation,Tdir
function [text] = readfilee(filepath)
%open file
fileID = fopen(filepath);
%read text by line
line=fgets(fileID); %read first line to ensure their is text
text=cell(1); %preallocate file input
i=1; %concatinates onto text cell so each line is its own cell

while line > -1 %while there is text on next line (last line is -1)
    if line(1) == '#' %if line starts with # do nothing
    else %save text input as text
        text{i}=line; %concatinate onto cell
        i=i+1;
        disp(line)
    end 
    line=fgets(fileID);
end