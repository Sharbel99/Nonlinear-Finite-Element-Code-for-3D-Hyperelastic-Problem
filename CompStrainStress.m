% Compute strain and stress at integration points
% Shoaib Goraya;10/2017

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

    %Extract patch nodal displacements
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
    
    [strain,stress] = CompStrainStress_Elem(xl,ul,mateprop,nel,ndf,elem);

    % Plot strain and stress arrays
%     fprintf('\nElement\t%1.0f\n\n',elem);
%     strain
%     stress
    
end