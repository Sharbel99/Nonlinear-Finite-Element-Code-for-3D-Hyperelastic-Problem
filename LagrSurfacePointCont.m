function [S_u_v,C_u_v] = LagrSurfacePointCont(P,C,nel,edge,r,s,dim)
%
% Copyright (C) Arif Masud and Tim Truster
%
%     -------------------------------------------------------------------

%   Compute Surface Point
    
    S_u_v = zeros(dim,1);
    C_u_v = 0;

    shp = shpl_2dl(r,s,nel);
    
    for dir = 1:dim
        for l = 1:nel
            S_u_v(dir) = S_u_v(dir) + shp(l)*P(l,dir);
        end
    end
    for l = 1:nel
        C_u_v = C_u_v + shp(l)*C(l);
    end
    
function shp = shpl_2dl(r,s,nel)

shp = zeros(1,nel);

if nel == 3
    
    shp(1,2)= r;
    shp(1,3)= s;
    shp(1,1)= 1-r-s;


elseif nel == 4

    shp(1,1) = 1/4*(1-r)*(1-s);
    shp(1,2) = 1/4*(1+r)*(1-s);
    shp(1,3) = 1/4*(1+r)*(1+s);
    shp(1,4) = 1/4*(1-r)*(1+s);
    
end