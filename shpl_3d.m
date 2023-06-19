function shp = shpl_3d(r,s,t,nel)
%
% Copyright (C) Arif Masud and Tim Truster
%
%
%           4 -------------- 3       2
%           |                |       | \
%           |                |       |  \
%           |                |       |   \
%           |                |       |    \
%           |                |       |     \
%           |                |       |      \
%           1 -------------- 2       3-------1
%

shp = zeros(4,nel);

if nel == 4
    
    shp(1,1)= 1;
    shp(2,1)= 0;
    shp(3,1)= 0;
    shp(4,1)= r;

    shp(1,2)= 0;
    shp(2,2)= 1;
    shp(3,2)= 0;
    shp(4,2)= s;

    shp(1,3)=-1;
    shp(2,3)=-1;
    shp(3,3)=-1;
    shp(4,3)= 1-r-s-t;

    shp(1,4)= 0;
    shp(2,4)= 0;
    shp(3,4)= 1;
    shp(4,4)= t;
    
elseif nel == 8
        ap1 = 1.0d0 + r;
        am1 = 1.0d0 - r;
        ap2 = 1.0d0 + s;
        am2 = 1.0d0 - s;
        ap3 = 1.0d0 + t;
        am3 = 1.0d0 - t;

%       Compute for ( - , - ) values

        c1      = 0.125*am1*am2;
        c2      = 0.125*am2*am3;
        c3      = 0.125*am1*am3;
        shp(1,1) = -c2;
        shp(1,2) =  c2;
        shp(2,1) = -c3;
        shp(2,4) =  c3;
        shp(3,1) = -c1;
        shp(3,5) =  c1;
        shp(4,1) =  c1*am3;
        shp(4,5) =  c1*ap3;

%       Compute for ( + , + ) values

        c1      = 0.125*ap1*ap2;
        c2      = 0.125*ap2*ap3;
        c3      = 0.125*ap1*ap3;
        shp(1,8) = -c2;
        shp(1,7) =  c2;
        shp(2,6) = -c3;
        shp(2,7) =  c3;
        shp(3,3) = -c1;
        shp(3,7) =  c1;
        shp(4,3) =  c1*am3;
        shp(4,7) =  c1*ap3;

%       Compute for ( - , + ) values

        c1      = 0.125*am1*ap2;
        c2      = 0.125*am2*ap3;
        c3      = 0.125*am1*ap3;
        shp(1,5) = -c2;
        shp(1,6) =  c2;
        shp(2,5) = -c3;
        shp(2,8) =  c3;
        shp(3,4) = -c1;
        shp(3,8) =  c1;
        shp(4,4) =  c1*am3;
        shp(4,8) =  c1*ap3;

%       Compute for ( + , - ) values

        c1      = 0.125*ap1*am2;
        c2      = 0.125*ap2*am3;
        c3      = 0.125*ap1*am3;
        shp(1,4) = -c2;
        shp(1,3) =  c2;
        shp(2,2) = -c3;
        shp(2,3) =  c3;
        shp(3,2) = -c1;
        shp(3,6) =  c1;
        shp(4,2) =  c1*am3;
        shp(4,6) =  c1*ap3 ;   
        
end    
end