% Assemble Element Stiffness into Model Stiffness
%
% Copyright (C) Arif Masud and Tim Truster
%
% 10/2009
% UIUC

locind1 = 0;
for ie1 = 1:nel
    for l1 = 1:ndf
        locind1 = locind1 + 1;
        grow = EDOFT(locind1);
        if grow <= neq
            locind2 = 0;
            for ie2 = 1:nel
                for l2 = 1:ndf
                    locind2 = locind2 + 1;
                    gcol = EDOFT(locind2);
                    if gcol <= neq
                        Kdd(grow, gcol) = Kdd(grow, gcol) + ElemK(locind1,locind2); %#ok<AGROW>
                    else %gcol > neq
                        gcol = gcol - neq;
                        Kdf(grow, gcol) = Kdf(grow, gcol) + ElemK(locind1,locind2); %#ok<AGROW>
                    end %if gcol
                end %l2
            end %ie2
        else %grow > neq
            grow = grow - neq;
            locind2 = 0;
            for ie2 = 1:nel
                for l2 = 1:ndf
                    locind2 = locind2 + 1;
                    gcol = EDOFT(locind2);
                    if gcol <= neq
                        Kfd(grow, gcol) = Kfd(grow, gcol) + ElemK(locind1,locind2); %#ok<AGROW>
                    else %gcol > neq
                        gcol = gcol - neq;
                        Kff(grow, gcol) = Kff(grow, gcol) + ElemK(locind1,locind2); %#ok<AGROW>
                    end %if gcol
                end %l2
            end %ie2
        end %if grow
    end %l1
end %ie1