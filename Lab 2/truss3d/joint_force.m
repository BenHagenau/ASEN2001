function [] = joint_force(connectivity,joints)

%calculate length of each member
%number of connections
x=size(connectivity);
%read connectivity line
for i=1:x(1)
    