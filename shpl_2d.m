function shp = shpl_2d(r,s,nel)
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

shp = zeros(3,nel);

if nel == 3
    
    shp(1,1)= 1;
    shp(2,1)= 0;
    shp(3,1)= r;

    shp(1,2)= 0;
    shp(2,2)= 1;
    shp(3,2)= s;

    shp(1,3)=-1;
    shp(2,3)=-1;
    shp(3,3)= 1-r-s;


elseif nel == 4
   
    % shp(i,j). j = denotes the node number,
    % i = 1, derivative w.r.t r
    % i = 2, derivative w.r.t s
    % i = 3, shape funtion
    
    shp(1,1) = -0.25*(1-s);
    shp(2,1) = -0.25*(1-r);
    shp(3,1) = 0.25*(1-r)*(1-s);
    
    shp(1,2) = 0.25*(1-s);
    shp(2,2) = -0.25*(1+r);
    shp(3,2) = 0.25*(1+r)*(1-s);
    
    shp(1,3) = 0.25*(1+s);
    shp(2,3) = 0.25*(1+r);
    shp(3,3) = 0.25*(1+r)*(1+s);
    
    shp(1,4) = -0.25*(1+s);
    shp(2,4) = 0.25*(1-r);
    shp(3,4) = 0.25*(1-r)*(1+s);
    
end