function [w,r,s] = intpntq(l,lint,ib)

% Copyright (C) Arif Masud and Tim Truster

%       integer l,lint,ib
%       real*8 ra(100),sa(100),wa(100),ls(9),lt(9),lw(9),w
%       real*8 x1d(10),w1d(10)
%       real*8 g4(4), h4(4)
%       integer l1,i,j,k
%       real*8 zero,pt1667,pt25,pt5,pt6,one,two,three,four,five,six 
%       real*8 etyone,g,h,forpt8,thirty,thrtsx,seven,r,s
% 
% 		                  
ls = [-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0,-1.d0,0.d0];
lt = [-1.d0,-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0,0.d0];
lw = [25.d0,25.d0,25.d0,25.d0,40.d0,40.d0,40.d0,40.d0,64.d0];
% c
%       data x1d /-0.973906528517172d0,-0.865063366688985d0,
%      &          -0.679409568299024d0,-0.433395394129247d0,
%      &          -0.148874338981631d0, 0.148874338981631d0,
%      &           0.433395394129247d0, 0.679409568299024d0,
%      &           0.865063366688985d0, 0.973906528517172d0/
% c
%       data w1d / 0.066671344308688d0, 0.149451349150581d0,
%      &           0.219086362515982d0, 0.269266719309996d0,
%      &           0.295524224714753d0, 0.295524224714753d0,
%      &           0.269266719309996d0, 0.219086362515982d0,
%      &           0.149451349150581d0, 0.066671344308688d0/
% 
% c.... define some constants (for feap)
        ra = zeros(100,1);
        sa = zeros(100,1);
        wa = zeros(100,1);
        zero   = 0.d0;
%         pt1667 = 0.1667d0
%         pt25   = 0.25d0
        pt5    = 0.5d0;
        one    = 1.0d0;
%         two    = 2.0d0;
        three  = 3.0d0;
        four   = 4.0d0;
        five   = 5.0d0;
        six    = 6.0d0;
% 
      if ib == 0
% 
      l1=sqrt(lint);
% c
% c     1x1 integration.
% c     ----------------
      if l1 == 1
      w     = four;
      ra(1) = zero;
      sa(1) = zero;
% !tjt      if (nel.eq.three) sa(1) = -one/three
      end    
% c
% c     2x2 integration.
% c     ----------------
      if l1 == 2
      g = one/sqrt(three);
      for i = 1:4
         w     = one;
         ra(i) = g*ls(i);
         sa(i) = g*lt(i);
      end
      end
% c
% c     3x3 integration.
% c     ----------------
      if l1 == 3
      pt6 = pt5+one/(two*five);
      etyone = three^four;
      g = sqrt(pt6);
      h = one/etyone;
      w = h*lw(l);
      for i = 1:9
         ra(i) = g*ls(i);
         sa(i) = g*lt(i);
      end
      end
% c
% c	4x4 integration.
% c	----------------
      if l1 == 4
	  forpt8 = four + four/five;
      thirty = five * six;
      thrtsx = six * six;
      seven  = one + six;
      g = sqrt(forpt8);
      h = sqrt(thirty)/thrtsx;
      g4(1) = sqrt((three + g)/seven);
      g4(4) = -g4(1);
      g4(2) = sqrt((three - g)/seven);
      g4(3) = -g4(2);
      h4(1) = pt5 - h;
      h4(2) = pt5 + h;
      h4(3) = pt5 + h;
      h4(4) = pt5 - h;
%
      i = zero;
      for j = 1:4
        for k = 1:4
	     i = i + 1;
	     wa(i)= h4(j)*h4(k);
	     ra(i)= g4(k);
	     sa(i)= g4(j);
        end
      end
       w = wa(l);
      end %if
% c
% c	10x10 integration.
% c	----------------
%       if (l1 .eq. 10) then
%       i = 1
%       do j = 1, 10
%       do k = 1, 10
% 	   ra(i) = x1d(k)
% 	   sa(i) = x1d(j)
% 	   wa(i) = w1d(k)*w1d(j)	   
% c	write(15,*)k,j,ra(i),sa(i),wa(i),i
%       i = i + 1
%       enddo
%       enddo
% c	rewind(15)
%       w = wa(l)
%       endif
% 
%       elseif (ib.eq.3) then
% c
% c	4x4 integration.
% c	----------------
%       if (lint .eq. 4) then
% 	  forpt8 = four + four/five
%       thirty = five * six
%       thrtsx = six * six
%       seven  = one + six
%       g = sqrt(forpt8)
%       h = sqrt(thirty)/thrtsx
%       g4(1) = sqrt((three + g)/seven)
%       g4(4) = -g4(1)
%       g4(2) = sqrt((three - g)/seven)
%       g4(3) = -g4(2)
%       h4(1) = pt5 - h
%       h4(2) = pt5 + h
%       h4(3) = pt5 + h
%       h4(4) = pt5 - h
% c
%       do 51 j = 1,4
% 	   wa(j)= h4(j)
% 	   ra(j)= -one
% 51	   sa(j)= g4(j)
%        w = wa(l)
%       endif
%       
%       elseif (ib.eq.4) then
% c
% c	4x4 integration.
% c	----------------
%       if (lint .eq. 4) then
% 	  forpt8 = four + four/five
%       thirty = five * six
%       thrtsx = six * six
%       seven  = one + six
%       g = sqrt(forpt8)
%       h = sqrt(thirty)/thrtsx
%       g4(1) = sqrt((three + g)/seven)
%       g4(4) = -g4(1)
%       g4(2) = sqrt((three - g)/seven)
%       g4(3) = -g4(2)
%       h4(1) = pt5 - h
%       h4(2) = pt5 + h
%       h4(3) = pt5 + h
%       h4(4) = pt5 - h
% c
%       do 61 j = 1,4
% 	   wa(j)= h4(j)
% 	   ra(j)= one
% 61	   sa(j)= g4(j)
%        w = wa(l)
%       endif
%       
%       elseif (ib.eq.1) then
% c
% c	4x4 integration.
% c	----------------
%       if (lint .eq. 4) then
% 	  forpt8 = four + four/five
%       thirty = five * six
%       thrtsx = six * six
%       seven  = one + six
%       g = sqrt(forpt8)
%       h = sqrt(thirty)/thrtsx
%       g4(1) = sqrt((three + g)/seven)
%       g4(4) = -g4(1)
%       g4(2) = sqrt((three - g)/seven)
%       g4(3) = -g4(2)
%       h4(1) = pt5 - h
%       h4(2) = pt5 + h
%       h4(3) = pt5 + h
%       h4(4) = pt5 - h
% c
%       do 71 j = 1,4
% 	   wa(j)= h4(j)
% 	   ra(j)= g4(j)
% 71	   sa(j)= -one
%        w = wa(l)
%       endif
%       
%       elseif (ib.eq.2) then
% c
% c	4x4 integration.
% c	----------------
%       if (lint .eq. 4) then
% 	  forpt8 = four + four/five
%       thirty = five * six
%       thrtsx = six * six
%       seven  = one + six
%       g = sqrt(forpt8)
%       h = sqrt(thirty)/thrtsx
%       g4(1) = sqrt((three + g)/seven)
%       g4(4) = -g4(1)
%       g4(2) = sqrt((three - g)/seven)
%       g4(3) = -g4(2)
%       h4(1) = pt5 - h
%       h4(2) = pt5 + h
%       h4(3) = pt5 + h
%       h4(4) = pt5 - h
% c
%       do 81 j = 1,4
% 	   wa(j)= h4(j)
% 	   ra(j)= g4(j)
% 81	   sa(j)= one
%        w = wa(l)
%       endif
%       
      end %if
% 
      r = ra(l);
      s = sa(l);
% 
%       return
% 
%       end