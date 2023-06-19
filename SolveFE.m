% Solve Partitioned Finite Element Matrix System
%
% Copyright (C) Arif Masud and Tim Truster
% 7/2009
% UIUC

if step==1 
    for i = 1:neq
        rhs = 0; 
        for j = 1:nieq 
            rhs = rhs + Kdf(i,j)*ModelDc(j)/NbLoadSteps; 
        end
        Fdtilda(i,1) = Fext(i,1) - rhs - Fint(i,1);
    end
else
    Fdtilda = Fext - Fint;
end

%Solve Kd = F

ModelDx = ModelDx+ Kdd\Fdtilda; %#ok<NASGU>