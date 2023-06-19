function ELDOFT = LocToGlobDOF(NFlags, NDOFT, maxnodes, ndir)
%
% Copyright (C) Arif Masud and Tim Truster
%
ELDOFT = zeros(1,ndir*maxnodes);
loc=0;

for node = 1:maxnodes
    
    for direction = 1:ndir
        loc = loc+1;
        ELDOFT(loc) = NDOFT(NFlags(node), direction);
    end
    
end