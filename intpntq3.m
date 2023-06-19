function [w,r,s,t] = intpntq(l,lint,ib)

% Copyright (C) Arif Masud and Tim Truster

% c
% c     1x1 integration.%%%%%%%%%%%%%Not updated
% c     ----------------
l1 = lint^(1/3);

      if l1 == 1
      wa(1) = 8;
      ra(1) = 0;
      sa(1) = 0;
      ta(1) = 0;
% !tjt      if (nel.eq.three) sa(1) = -one/three
      end    
% c
% c     2x2x2 integration.
% c     ----------------
      if l1 == 2
      g(1) = -0.577350269189626;
      g(2) = 0.577350269189626;
      n=0;
      for i = 1:2
          for j=1:2
              for k=1:2
                  n=n+1;
                  ra(n) = g(k);
                  sa(n) = g(j);
                  ta(n) = g(i);
                  wa(n) = 1;
              end
          end
      end
      
      end
% c    %%%%%%%%%%%%  everything from here onwards is unchanged
% c     3x3 integration.
% c     ----------------
      if l1 == 3
      g(1)= -0.774596669241483;
      g(2)= 0;
      g(3)= 0.774596669241483;
      weight(1)= 0.555555555555556;
      weight(2)= 0.888888888888889;
      weight(3)= 0.555555555555556;
      n=0;
      for i=1:3
          for j=1:3
              for k=1:3
                  n=n+1;
                  ra(n)=g(k);
                  sa(n)=g(j);
                  ta(n)=g(i);
                  wa(n)=weight(i)*weight(j)*weight(k);
              end
          end
      end
      end
% c
% c	4x4 integration.
% c	----------------
      if l1 == 4
      g(1)= -0.861136311594953;
      g(2)= -0.339981043584856;
      g(3)= 0.339981043584856;
      g(4)= 0.861136311594953;
      weight(1)= 0.347854845137454;
      weight(2)= 0.652145154862546;
      weight(3)= 0.652145154862546;
      weight(4)= 0.347854845137454;
      n=0;
      for i=1:4
          for j=1:4
              for k=1:4
                  n=n+1;
                  ra(n)=g(k);
                  sa(n)=g(j);
                  ta(n)=g(i);
                  wa(n)=weight(i)*weight(j)*weight(k);
              end
          end
      end
      end %if

% 
      r = ra(l);
      s = sa(l);
      t = ta(l);
      w = wa(l);
end