function [cartd1, xsj] = shpg_2d(shp,xl,nel)
%      
%      * * F E A P * * A Finite Element Analysis Program
%
%....  Copyright (c) 1984-2002: Regents of the University of California
%                               All rights reserved
%
%....  Modified for NURBS functions by Tim Truster, UIUC,  05/25/2009
%....  Corrected to match ordering of Hughes and Masud,    07/11/2009
%
%-----[--.----+----.----+----.-----------------------------------------]
%      Purpose:
%
%      Inputs:
%         ss(*)    - Gauss points
%
%      Outputs:
%         shp(4,*) - Quadratic shape functions
%-----[--.----+----.----+----.-----------------------------------------]

%      integer    i,nel2,ir(27),is(27),it(27),ndm,j,k,inode,
%     &		ideriv,n

%      real*8     nr(3),dr(3),ns(3),ds(3),nt(3),dt(3),ss(3),shp(4,nel2+1)
%     &		,c1(6,3),tc(6,3),tcj(6,3),at1(6),at2(6),xsj,d1,d2,d3,rxsj,
%     &		t2(6,6),cartd1(4,nel2+1),cartd2(6,nel2+1),shp2(6,nel2+1),
%     &		xl(ndm,nel2-1),xs(3,3),ad(3,3),x,y,z,ddr(3),dds(3),ddt(3)

xs = zeros(2,2);
rxs = zeros(2,2);
cartd1 = zeros(3,nel);
ad = zeros(2,2);
		  
% compute jacobian transformation

% first derivatives only

        for deriv = 1:2
            for xyz = 1:2
%               xs(i,j)=0;
                for inode = 1:nel
                    xs(xyz,deriv)=xs(xyz,deriv)+shp(deriv,inode)*xl(xyz,inode);
                end % inode
            end % xyz
        end %deriv
%         xs = shp(1:2,1:nel)*xl(1:2,1:nel)';

%     Compute adjoint to jacobian

      ad(1,1) = xs(2,2);
      ad(1,2) = -xs(1,2);

      ad(2,1) = -xs(2,1);
      ad(2,2) = xs(1,1);

%     Compute determinant of jacobian

      xsj  = xs(1,1)*xs(2,2)- xs(1,2)*xs(2,1);
%       xsj = xs(1,1:2)*ad(1:2,1);

      rxsj = 1/xsj;
      
%     Compute jacobian inverse

      rxs(1,1) = ad(1,1)*rxsj;
      rxs(2,2) = ad(2,2)*rxsj;
      rxs(1,2) = ad(1,2)*rxsj;
      rxs(2,1) = ad(2,1)*rxsj;

%     Compute derivatives with repect to global coords.

      for k = 1:nel
          d1 = shp(1,k)*rxs(1,1) + shp(2,k)*rxs(2,1);
          d2 = shp(1,k)*rxs(1,2) + shp(2,k)*rxs(2,2);
          cartd1(1,k) = d1;
          cartd1(2,k) = d2;
          cartd1(3,k) = shp(3,k);
      end
%       cartd1 = [rxs(1:2,1:2)'*shp(1:2,1:nel)
%                 shp(3,1:nel)];
