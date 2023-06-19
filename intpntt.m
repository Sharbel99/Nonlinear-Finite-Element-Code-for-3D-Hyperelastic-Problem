function [w,c1,c2] = intpntt(l,lint,ib)
% 
% Copyright (C) Arif Masud and Tim Truster
%
% 	implicit none
% c.... remove above card for single precision operation
% c
% 	integer l,lint,ib
% 	integer i
% 	real*8 zero,pt1667,pt25,pt5,one,two,three,four,five,six 
% 	real*8 w,w1,r1,w2,r2,w3a,w3b,r3a,r3b1,r3b2,w7a,w7b,w7d,r7a,r7c,
%      &	   r7b,r7e,r7d,w8a,w8c,w8e,r8b,r8a,r8d,r8c,r8e,r8f,r8g
% 	real w9a,w9b,w9c,w9d,r9a,r9b,r9c,r9d,r9e,r9f,r9g,r9h
% 	real*8 forpt8,thirty,thrtsx,seven,g4(4),h4(4),g,h
%       real*8 cl1(16),cl2(16),cl3(16),c1,c2
% 	integer xsi1,xsi2,xsi3
% 
% c

cl1 = zeros(16,1);
cl2 = zeros(16,1);
cl3 = zeros(16,1);

r1 = 0.33333333333333333333d00;
w1 = 1.d00;
% r2 = 0.5d00;
% w2 = 0.3333333333333333333d00;
r3a = 0.3333333333333333333d00;
w3a = -0.5625d00;
r3b1 = 0.6d00;
w3b = 0.520833333333333333d00;
r3b2 = 0.2d00;
r7a = 0.3333333333333333333d00;
w7a = 0.2250300003d00;
r7b = 0.05971587178977d00;
w7b = 0.132394152788506d00;
r7c = 0.470142064105115d00;
r7d = 0.797426985353087d00;
w7d = 0.125939180544827d00;
r7e = 0.101286507323456d00;
% r8a = 0.873821971016996d00;
% w8a = 0.050844906370207d00;
% r8b = 0.063089014491502d00;
% r8c = 0.501426509658179d00;
% w8c = 0.116786275726379d00;
% r8d = 0.249286745170910d00;
% r8e = 0.636502499121399d00;
% w8e = 0.082851075618374d00;
% r8f = 0.310352451033785d00;
% r8g = 0.053145049844816d00;
% c
r9a = 0.0651301029022d00;
w9a = 0.0533472356088d00;
r9b = 0.8697397941956d00;
r9c = 0.3128654960049d00;
w9b = 0.0771137608903d00;
r9d = 0.6384441885698d00;
r9e = 0.0486903154253d00;
r9f = 0.2603459660790d00;
w9c = 0.1756152574332d00;
r9g = 0.4793080678419d00;
w9d = -0.1495700444677d00;
r9h = 0.3333333333333d00;
 
% 	   common /jumpdataT/ xsi1,xsi2,xsi3
% c
% c.... define some constants (for feap)
zero   = 0.d0;
% pt1667 = 0.1667d0;
% pt25   = 0.25d0;
pt5    = 0.5d0;
one    = 1.0d0;
two    = 2.0d0;
three  = 3.0d0;
four   = 4.0d0;
five   = 5.0d0;
six    = 6.0d0;

if ib == 0
        
% c
% c     One point integration
% c     ---------------------
    if lint == 1
        w=w1/two;
        cl1(1)=r1;
        cl2(1)=r1;
        cl3(1)=one-r1-r1;
    end
% c
% c     Three point integration
% c     -----------------------
%         if lint == 3
% c            
%             w=w2/two;
% 
% 
% 	  if(xsi1.eq.0) then
% 
% 		  cl1(1)=zero					!r1=0
%             cl2(1)=r2						!s1=0.5
%             cl3(1)=one-cl1(1)-cl2(1)		!t1=1-r1-s1
% 
%             cl1(2)=zero					!r2=0
%             cl2(2)=r2						!s2=0.5
%             cl3(2)=one-cl1(2)-cl2(2)		!t2=1-r2-s2
% 
%             cl1(3)=zero					!r3=0
%             cl2(3)=zero					!s3=0.0
%             cl3(3)=one-cl1(3)-cl2(3)		!t3=1-r3-s3
% 
% 	   
% 
% 	  elseif(xsi2.eq.0) then
% 
% 
% 		  cl1(1)=r2						!r1=0.5
%             cl2(1)=zero					!s1=0
%             cl3(1)=one-cl1(1)-cl2(1)		!t1=1-r1-s1
% 
%             cl1(2)=zero					!r2=0.0
%             cl2(2)=zero					!s2=0
%             cl3(2)=one-cl1(2)-cl2(2)		!t2=1-r2-s2
% 
%             cl1(3)=r2						!r3=0.5
%             cl2(3)=zero					!s3=0
%             cl3(3)=one-cl1(3)-cl2(3)		!t3=1-r3-s3
% 
% 
% 
% 	  elseif(xsi3.eq.0) then
% 
% 
% 		  cl1(1)=r2						!r1=0.5
%             cl2(1)=one-cl1(1)				!s1=1-r1
%             cl3(1)=zero					!t1=0
% 
%             cl1(2)=zero					!r2=0.0
%             cl2(2)=one-cl1(2)				!s2=1-r2
%             cl3(2)=zero					!t2=0
% 
%             cl1(3)=r2						!r3=0.5
%             cl2(3)=one-cl1(3)				!s3=1-r3
%             cl3(3)=zero					!t3=0
% 
% 
% 	  endif	! On xsi1,xsi2,xsi3
% 
% 
%         end if  ! on lint.eq.3
% 
% 
% cjp	if(lint.eq.3) then
% c            
% cjp            w=w2/two
% cjp            cl1(1)=r2		
% cjp            cl2(1)=r2		
% cjp            cl3(1)=zero
% 	
% cjp            cl1(2)=zero
% cjp            cl2(2)=r2
% cjp            cl3(2)=r2
% 
% cjp            cl1(3)=r2
% cjp            cl2(3)=zero
% cjp            cl3(3)=r2
% cjp        end if

% c
% c     Four point integration
% c     ----------------------
    if lint == 4
        if l == 1 
            w= w3a/two;
        end
        if l>=2 && l<=4
            w=w3b/two;
        end
        cl1(1)=r3a;
        cl2(1)=r3a;

% cjp	cl3(1)=1.d00 - r3a - r3a  T6 sqrt(-) with lint=4

        cl3(1)=one - r3a - r3a;
        cl1(2)=r3b1;
        cl2(2)=r3b2;
        cl3(2)=r3b2;
        cl1(3)=r3b2;
        cl2(3)=r3b1;
        cl3(3)=r3b2;
        cl1(4)=r3b2;
        cl2(4)=r3b2;
        cl3(4)=r3b1;
    end
% c
% c     Seven point integration
% c     -----------------------
    if lint == 7 
        if l == 1  
            w= w7a/two;
        end
        if l>=2 && l<=4
            w=w7b/two;
        end
        if l>=5 && l<=7
            w=w7d/two;
        end
        cl1(1)=r7a;
        cl2(1)=r7a;
        cl3(1)=r7a;
        for i=2:4
            cl1(i)=r7c;
            cl2(i)=r7c;
            cl3(i)=r7c;
        end
        cl1(2)=r7b;
        cl2(3)=r7b;
        cl3(4)=r7b;
        for i=5:7
            cl1(i)=r7e;
            cl2(i)=r7e;
            cl3(i)=r7e;
        end 
        cl1(5)=r7d;
        cl2(6)=r7d;
        cl3(7)=r7d;
    end
% c 
% c	12 point integration
% c	--------------------
% 	  if lint == 12
% 
%           if (l>=1 && l<=3)  
%               w=w8a/two;
%           end
%           if (l>=4 && l<=6)
%               w=w8c/two;
%           end
%           if (l>=7 && l<=12)
%               w=w8e/two;
%           end
% 
% 
% 	  if(xsi1.eq.0) then
% 
% 	  do i=1,12
% 
% 	   cl1(i)=0.d0
% 
%       enddo
% 
% 		cl2(1)=r8b
% 		cl3(1)=1-cl1(1)-cl2(1)
% 
% 		cl2(2)=r8a
% 		cl3(2)=1-cl1(2)-cl2(2)
% 
% 		cl2(3)=r8b
% 		cl3(3)=1-cl1(3)-cl2(3)
% 	
% 		cl2(4)=r8d
% 		cl3(4)=1-cl1(4)-cl2(4)
% 
% 		cl2(5)=r8c
% 		cl3(5)=1-cl1(5)-cl2(5)
% 
% 		cl2(6)=r8d
% 		cl3(6)=1-cl1(6)-cl2(6)
% 
% 		cl2(7)=r8f
% 		cl3(7)=1-cl1(7)-cl2(7)
% 
% 		cl2(8)=r8e
% 		cl3(8)=1-cl1(8)-cl2(8)
% 
%           cl2(9)=r8e
%           cl3(9)=1-cl1(9)-cl2(9)
% 
%           cl2(10)=r8f
%           cl3(10)=1-cl1(10)-cl2(10)
% 
%           cl2(11)=r8g
%           cl3(11)=1-cl1(11)-cl2(11)
% 
%           cl2(12)=r8g
%           cl3(12)=1-cl1(12)-cl2(12)
% 
% 	
% 
% 	  elseif(xsi2.eq.0) then
% 
% 	  do i=1,12
% 
% 	   cl2(i)=0.d0
% 
%       enddo
% 
% 
% 		cl1(1)=r8a
% 		cl3(1)=1-cl1(1)-cl2(1)
% c
% 		cl1(2)=r8b
% 		cl3(2)=1-cl1(2)-cl2(2)
% c
% 		cl1(3)=r8b
% 		cl3(3)=1-cl1(3)-cl2(3)
% c		
% 		cl1(4)=r8c
% 		cl3(4)=1-cl1(4)-cl2(4)
% c
% 		cl1(5)=r8d
% 		cl3(5)=1-cl1(5)-cl2(5)
% c
% 		cl1(6)=r8d
% 		cl3(6)=1-cl1(6)-cl2(6)
% c
% 		cl1(7)=r8e
% 		cl3(7)=1-cl1(7)-cl2(7)
% c
% 		cl1(8)=r8f
% 		cl3(8)=1-cl1(8)-cl2(8)
% c
%           cl1(9)=r8g
%           cl3(9)=1-cl1(9)-cl2(9)
% c
%           cl1(10)=r8g
%           cl3(10)=1-cl1(10)-cl2(10)
% c
%           cl1(11)=r8f
%           cl3(11)=1-cl1(11)-cl2(11)
% c
%           cl1(12)=r8e
%           cl3(12)=1-cl1(12)-cl2(12)
% c
% 
% 
% 	  elseif(xsi3.eq.0) then
% 
% 	  do i=1,12
% 
% 	   cl3(i)=0.d0
% 
%       enddo
% 
% 		cl1(1)=r8a
% 		cl2(1)=1-cl1(1)
% c
% 		cl1(2)=r8b
% 		cl2(2)=1-cl1(2)
% c
% 		cl1(3)=r8b
% 		cl2(3)=1-cl1(3)
% c		
% 		cl1(4)=r8c
% 		cl2(4)=1-cl1(4)
% c
% 		cl1(5)=r8d
% 		cl2(5)=1-cl1(5)
% c
% 		cl1(6)=r8d
% 		cl2(6)=1-cl1(6)
% c
% 		cl1(7)=r8e
% 		cl2(7)=1-cl1(7)
% c
% 		cl1(8)=r8f
% 		cl2(8)=1-cl1(8)
% c
%           cl1(9)=r8g
%           cl2(9)=1-cl1(9)
% c
%           cl1(10)=r8g
%           cl2(10)=1-cl1(10)
% c
%           cl1(11)=r8f
%           cl2(11)=1-cl1(11)
% c
%           cl1(12)=r8e
%           cl2(12)=1-cl1(12)
% c
% 
% 
% 	  endif	! On xsi1,xsi2,xsi3
% 
%      
%       end if	! On lint=12


% c
% c	13 point integration
% c	--------------------
    if lint == 13

        if (l >= 1 && l <= 3)
          w=w9a/two;
        end
        if (l >= 4 && l <= 9)
          w=w9b/two;
        end
        if (l >= 10 && l <= 12)
          w=w9c/two;
        end
        if (l == 13)
          w=w9d/two;
        end

        cl1(1)=r9a;
        cl1(2)=r9b;
        cl1(3)=r9a;
        cl1(4)=r9c;
        cl1(5)=r9d;
        cl1(6)=r9e;
        cl1(7)=r9d;
        cl1(8)=r9c;
        cl1(9)=r9e;
        cl1(10)=r9f;
        cl1(11)=r9g;
        cl1(12)=r9f;
        cl1(13)=r9h;

        cl2(1)=cl1(1);
        cl2(2)=cl1(1);
        cl2(3)=cl1(2);
        cl2(4)=cl1(6);
        cl2(5)=cl1(4);
        cl2(6)=cl1(5);
        cl2(7)=cl1(6);
        cl2(8)=cl1(5);
        cl2(9)=cl1(4);
        cl2(10)=cl1(10);
        cl2(11)=cl1(10);
        cl2(12)=cl1(11);
        cl2(13)=cl1(13);

        for i=1:13
          cl3(i)=1.d0-cl1(i)-cl2(i);
        end

    end

elseif ib == 1
% c
% c	4x4 integration.
% c	----------------
    if (lint == 4)
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

        cl1(l)= (g4(l)+one)/two;
        cl2(l)= zero;
        w = h4(l)/2; %triangle coords are half of Guass coords
    end

% elseif (ib.eq.2) then
% 
% elseif (ib.eq.3) then

end %ib

c1 = cl1(l);
c2 = cl2(l);


end