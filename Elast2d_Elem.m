function [ElemK, ElemF] = Elast2d_Elem(xl,mateprop,nel,ndf,ndm,ul)
%
% Copyright (C) Arif Masud and Tim Truster
%
% Subroutine to compute stiffness matrix and force vector for linear
% 2-dimensional elasticity element. Element currently supports bilinear
% quadrilateral elements with the following node and shape function
% labelling scheme:
%
%  (-1, 1)  4 -------------- 3 ( 1, 1)
%           |       s        |
%           |       ^        |
%           |       |        |
%           |       .-> r    |
%           |                |
%           |                |
%  (-1,-1)  1 -------------- 2 ( 1,-1)
%
% Element local coordinates (r,s) are defined by a coordinate axis with the
% origin at the center of the element; the corners of the element have
% local coordinate values as shown in the figure.
%
% Definitions for input:
%
%   xl:              = local array containing (x,y) coordinates of nodes
%                      forming the element; format is as follows:
%                          Nodes    |        n1  n2  n3  n4
%                          x-coord  |  xl = [x1  x2  x3  x4
%                          y-coord  |        y1  y2  y3  y4];
%
%   mateprop:        = vector of material properties:
%                          mateprop = [E v t]; 
%                                   = [(Young's Modulus) (Poisson's Ratio)
%                                      (thickness)];
%
%   nel:             = number of nodes on current element (4)
%
%   ndf:             = max number of DOF per node (2)
%
%   ndm:             = space dimension of mesh (2)
%
%   ul:            = displacement of nodes
%
% Definitions for output:
%
%   ElemK:           = element stiffness matrix containing stiffness
%                      entries in the following arrangement, where
%                      wij corresponds to weighting function (i), coordinate
%                      direction (j), and ukl corresponds to displacement
%                      function (k), coordinate direction (l):
%                                 u1x  u1y  u2x  u2y  u3x  u3y  u4x  u4y
%                      w1x  ElemK[ .    .    .    .    .    .    .    .
%                      w1y         .    .    .    .    .    .    .    .
%                      w2x         .    .    .    .    .    .    .    .
%                      w2y         .    .    .    .    .    .    .    .
%                      w3x         .    .    .    .    .    .    .    .
%                      w3y         .    .    .    .    .    .    .    .
%                      w4x         .    .    .    .    .    .    .    .
%                      w4y         .    .    .    .    .    .    .    . ];
%
%   ElemF:           = element force vector containing force entries in the
%                      following arrangement:
%                      w1x  ElemF[ . 
%                      w1y         . 
%                      w2x         . 
%                      w2y         . 
%                      w3x         . 
%                      w3y         . 
%                      w4x         . 
%                      w4y         . ];                      
%
% Definitions of local constants:
%
%   nst:             = size of element arrays (ndf*nel)
%
%



% Set Material Properties

Elema = mateprop(1);
Elemb = mateprop(2);
thick = mateprop(3);
ul_elem = reshape(ul,ndf*nel,1);

% Initialize Matrix and Vector

nst = nel*ndf;
ElemK = zeros(nst);
ElemF = zeros(nst,1);

% Load Guass Integration Points

if nel == 3
    lint = 4;
else
    lint = 4;
end

% Loop over integration points
for l = 1:lint

        if nel == 3
            [Wgt,r,s] =  intpntt(l,lint,0);
        else
            [Wgt,r,s] =  intpntq(l,lint,0);
        end

        % Evaluate local basis functions at integration point
        shp = shpl_2d(r,s,nel);

        % Evaluate first derivatives of basis functions at int. point
        [Qxy, Jdet] = shpg_2d(shp,xl,nel);

        % Form B matrix
        if nel == 3
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3)];
        else
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        Qxy(1,4) 0 
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3) 0        Qxy(2,4)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3) Qxy(2,4) Qxy(1,4)];  
        end


        
        % Material moduli contribution of D
        
        Dmat=Elema./(1-Elemb.^2).*[1         Elemb   0
                                                 Elemb  1         0        
                                                 0        0         (1-Elemb)/2];
           

        
        % Update integration weighting factor
        W = Wgt*Jdet*thick;

        ElemK = ElemK + W*Bmat'*Dmat*Bmat;

end %je