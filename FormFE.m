% Assemble Quantities from Model Routine
%
% Copyright (C) Arif Masud and Tim Truster
%
% 7/2009
% UIUC

Kdd = zeros(neq,neq);
Kdf = zeros(neq,nieq);
Kfd = zeros(nieq,neq);
Kff = zeros(nieq,nieq);

AssemQuant = 'AssemStifForc';
% AssemQuant = 'AssemStifForc';

for elem = 1:numel
    
    %Determine element size parameters
  if ndm ==2
    if nen == 3
        nel = 3;
    elseif nen == 4
        if ix(elem,nen) == 0
            nel = 3;
        else
            nel = 4;
        end
    elseif nen == 6
        nel = 6;
    else
        if ix(elem,nen) == 0
            nel = 6;
        else
            nel = 9;
        end
    end
  else 
    if nen == 4
        nel = 4;
    elseif nen == 8
        if ix(elem,nen) == 0
            nel = 4;
        else
            nel = 8;
        end
    end
  end
    nst = nel*ndf;
    
    %Extract patch nodal coordinates
    xl = zeros(ndm, nel);
    ElemFlag = zeros(nel, 1);
    for k = 1:nel
        node = ix(elem,k);
        ElemFlag(k) = node;
        for l = 1:ndm
            xl(l,k) = NodeTable(node,l);
        end
    end
    
    %Extract patch material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);
    
    %Compute and Assemble Patch Stiffness
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    
    ul = zeros(ndm,nel);
    for i = 1:nel*ndf
        ndof_index = EDOFT(i);
        if(ndof_index<=neq)
            ul(i) = ModelDx(ndof_index);
        else
            ul(i) = ModelDc(ndof_index-neq);
        end
    end
    
    switch iel
        case 1 %Small-Deformation Isotropic Elastostatics Element
            if ndm == 3
                [ElemK, ElemF,C_ep(:,:,:,step+1),stress,strain] = ...
                    Elast3d_Elem(xl,mateprop,nel,ndf,ul);
            else %ndm == 2
                [ElemK, ElemF] = Elast2d_Elem(xl,mateprop,nel,ndf,ndm,ul);
            end
        case 2
            if ndm == 3

            else %ndm == 2
                L_Elem2_2d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% should add mine, very similar to Elst2d
            end
        case 3 %Stabilized Mixed Pressure-Displacement Element
            if ndm == 2
                L_Elem3_2d
            else %ndm == 3

            end
        case 4 %Implicit Error Element
            if ndm == 2
                L_Elem4_2d
            else %ndm == 3

            end
        case 5 %Stabilized Mixed Pressure-Displacement Element, Error
            if ndm == 2
                L_Elem5_2d
            else %ndm == 3

            end
    end

    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end