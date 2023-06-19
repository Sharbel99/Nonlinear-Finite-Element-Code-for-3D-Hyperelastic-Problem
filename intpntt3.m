function [w,c1,c2,c3] = intpntt3(l,lint,ib)


if ib == 0
        
% c
% c     One point integration
% c     ---------------------
    if lint == 1
        w1(1)  = 1/6;
        cl1(1) = 0.25;
        cl2(1) = 0.25;
        cl3(1) = 0.25;
    end
% c

% c     Four points integration
% c     ---------------------
    if lint == 4
        for i=1:4
            cl1(i) = 0.1381966011250105;
            cl2(i) = 0.1381966011250105;
            cl3(i) = 0.1381966011250105;
            w1(i)  = 0.25/6;
        end
        cl1(1) = 0.5854101966249658;
        cl2(2) = cl1(1);
        cl3(i) = cl1(1);
    end
% c

end %ib

c1 = cl1(l);
c2 = cl2(l);
c3 = cl3(l);
w  = w1(l);

end