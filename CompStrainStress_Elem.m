function [strain,stress] = CompStrainStress_Elem(xl,ul,mateprop,nel,ndf,elem)
%
% Subroutine to compute strain and stress for linear
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
%   PSPS:            = flag for plane stress ('s') or plane strain ('n')
%
% Definitions for output:
%
%   strain:          = strain array containing strain components
%                      at integration points:
%                                    xx   yy   xy
%                      int1  strain[ .    .    .   
%                      int2          .    .    .  
%                      int3          .    .    .  
%                      int4          .    .    .  ];
%
%   stress:          = stress array containing stress components
%                      at integration points:
%                                    xx   yy   xy
%                      int1  stress[ .    .    .   
%                      int2          .    .    .  
%                      int3          .    .    .  
%                      int4          .    .    .  ];                    
%
% Definitions of local constants:
%
%   nst:             = size of element arrays (ndf*nel)
%
%

% Set Material Properties

% Set Material Properties

ElemE = mateprop(1);
Elemv = mateprop(2);


% Load Guass Integration Points

if nel == 4
    lint = 4; 
else
    lint = 8; 
end

A=ElemE*(1-Elemv)/((1+Elemv)*(1-2*Elemv));
B=ElemE*Elemv/((1+Elemv)*(1-2*Elemv));
G=ElemE/2/(1+Elemv);

    
Dmat = [A B B 0 0 0
            B A B 0 0 0
            B B A 0 0 0
            0 0 0 G 0 0
            0 0 0 0 G 0
            0 0 0 0 0 G];
    
% Initialize Matrix and Vector

nst = nel*ndf;
strain = zeros(lint,6);
stress = zeros(lint,6);
ul_elem = reshape(ul,ndf*nel,1);

% Loop over integration points
for l = 1:lint

        if nel == 4
            [Wgt,r,s,t] =  intpntt3(l,lint,0);
        else
            [Wgt,r,s,t] =  intpntq3(l,lint,0);
        end

        % Evaluate local basis functions at integration point
        shp = shpl_3d(r,s,t,nel);

        % Evaluate first derivatives of basis functions at int. point
        [Qxy, Jdet] = shpg_3d(shp,xl,nel);

        % Form B matrix
        if nel == 4
        Bmat = [Qxy(1,1) 0        0        Qxy(1,2) 0        0       ,...
                    Qxy(1,3) 0        0        Qxy(1,4) 0        0       ;...
                    0        Qxy(2,1) 0        0        Qxy(2,2) 0       ,...
                    0        Qxy(2,3) 0        0        Qxy(2,4) 0       ;...
                    0        0        Qxy(3,1) 0        0        Qxy(3,2),...
                    0        0        Qxy(3,3) 0        0        Qxy(3,4);...
                    Qxy(2,1) Qxy(1,1) 0        Qxy(2,2) Qxy(1,2) 0       ,...
                    Qxy(2,3) Qxy(1,3) 0        Qxy(2,4) Qxy(1,4) 0       ;...
                    0        Qxy(3,1) Qxy(2,1) 0        Qxy(3,2) Qxy(2,2),...
                    0        Qxy(3,3) Qxy(2,3) 0        Qxy(3,4) Qxy(2,4);...
                    Qxy(3,1) 0        Qxy(1,1) Qxy(3,2) 0        Qxy(1,2),...
                    Qxy(3,3) 0        Qxy(1,3) Qxy(3,4) 0        Qxy(1,4)];
        else
        Bmat = [Qxy(1,1) 0        0        Qxy(1,2) 0        0       ,...
                    Qxy(1,3) 0        0        Qxy(1,4) 0        0       ,...
                    Qxy(1,5) 0        0        Qxy(1,6) 0        0       ,...
                    Qxy(1,7) 0        0        Qxy(1,8) 0        0       ;...
                    0        Qxy(2,1) 0        0        Qxy(2,2) 0       ,...
                    0        Qxy(2,3) 0        0        Qxy(2,4) 0       ,...
                    0        Qxy(2,5) 0        0        Qxy(2,6) 0       ,...
                    0        Qxy(2,7) 0        0        Qxy(2,8) 0       ;...
                    0        0        Qxy(3,1) 0        0        Qxy(3,2),...
                    0        0        Qxy(3,3) 0        0        Qxy(3,4),...
                    0        0        Qxy(3,5) 0        0        Qxy(3,6),...
                    0        0        Qxy(3,7) 0        0        Qxy(3,8);...
                    Qxy(2,1) Qxy(1,1) 0        Qxy(2,2) Qxy(1,2) 0       ,...
                    Qxy(2,3) Qxy(1,3) 0        Qxy(2,4) Qxy(1,4) 0       ,...
                    Qxy(2,5) Qxy(1,5) 0        Qxy(2,6) Qxy(1,6) 0       ,...
                    Qxy(2,7) Qxy(1,7) 0        Qxy(2,8) Qxy(1,8) 0       ;...
                    0        Qxy(3,1) Qxy(2,1) 0        Qxy(3,2) Qxy(2,2),...
                    0        Qxy(3,3) Qxy(2,3) 0        Qxy(3,4) Qxy(2,4),...
                    0        Qxy(3,5) Qxy(2,5) 0        Qxy(3,6) Qxy(2,6),...
                    0        Qxy(3,7) Qxy(2,7) 0        Qxy(3,8) Qxy(2,8);...
                    Qxy(3,1) 0        Qxy(1,1) Qxy(3,2) 0        Qxy(1,2),...
                    Qxy(3,3) 0        Qxy(1,3) Qxy(3,4) 0        Qxy(1,4),...
                    Qxy(3,5) 0        Qxy(1,5) Qxy(3,6) 0        Qxy(1,6),...
                    Qxy(3,7) 0        Qxy(1,7) Qxy(3,8) 0        Qxy(1,8)];           
        end
        
        % get strain
        strain_temp = Bmat*ul_elem;
        
        strain(l,1) = strain_temp(1);
        strain(l,2) = strain_temp(2);
        strain(l,3) = strain_temp(3);
        strain(l,4) = strain_temp(4)/2;
        strain(l,5) = strain_temp(5)/2;
        strain(l,6) = strain_temp(6)/2;
           
        % get stress
        stress_temp=Dmat*strain_temp;
        for i=1:6
            stress(l,i)=stress_temp(i);
        end

end