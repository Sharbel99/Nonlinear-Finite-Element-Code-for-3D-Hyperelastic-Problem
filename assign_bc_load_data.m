% function [ModelFc,ModelDc,neq,nieq,NDOFT] = assign_bc_load_data(numnp, NodeTable, PatchTable, PatchIndex, FaceBC, NodeBC, FaceLoad, NodeLoad, BCLIndex)
%
% Copyright (C) Arif Masud and Tim Truster
%
% 03/28/2009

NDOFT = zeros(numnp,2*ndf);
nieq = 0;

%Get BC
len = BCLIndex(1);
for i = 1:len
    node = NodeBC(i,1);
    dir = NodeBC(i,2);
    displacement = NodeBC(i,3);
    NDOFT(node,dir) = -1;
    NDOFT(node,dir+ndf) = displacement;
    nieq = nieq + 1;
end
%Assign DOF
neq = ndf*numnp - nieq;
ModelDc = zeros(nieq, 1);
na = 0;
ni = neq; 
for i = 1:numnp
    for j = 1:ndf
        if NDOFT(i,j) == 0
            na = na + 1;
            NDOFT(i,j) = na;
        elseif NDOFT(i,j) == -1
            ni = ni + 1;
            NDOFT(i,j) = ni;
            ModelDc(ni-neq) = NDOFT(i,j+ndf);
        else
            NDOFT(i,j) = NDOFT(NDOFT(i,j),j);
        end
    end
end

%Get Loads
Fd = zeros(neq, 1);

len = BCLIndex(2);
for i = 1:len
    node = NodeLoad(i,1);
    dir = NodeLoad(i,2);
    force = NodeLoad(i,3);
    dof = NDOFT(node,dir);
    if dof <= neq
        Fd(dof) = force;
    end
end