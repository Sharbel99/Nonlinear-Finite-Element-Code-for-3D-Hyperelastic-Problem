% Linear Finite Element Program for 2-D Elasticity
%
% Copyright (C) Arif Masud and Tim Truster
%
% This program computes a numerical solution to a finite element model
% using input on the geometry and physical properties of a mesh, and on the
% boundary conditions and applied loads. The routine assembles element
% quantities into the stiffness matrix and force vector to create a system
% of equations which is then solved for the nodal values of the
% displacement field. Boundary conditions are applied to constrain the
% stiffness matrix and to augment the force vector. The output is a list of
% the displacements printed on screen, contour plots of the displacement
% fields, and a plot of the deformed configuration of the mesh.
%
% Mesh input should be uploaded by running an input .m file before
% executing this program.
%
% Format of required input:
%
%   numnp:           = number of nodes in the mesh (length(NodeTable))
%
%   numel:           = number of elements in the mesh
%
%   nen:             = maximum number of nodes per element (4)
%
%   PSPS:            = flag for plane stress ('s') or plane strain ('n')
%
%   NodeTable:       = table of mesh nodal coordinates defining the
%                      geometry of the mesh; format of the table is as
%                      follows:
%                          Nodes  |             x-coord  y-coord
%                          n1     |  NodeTable = [x1     y1
%                          n2     |               x2     y2
%                          ...    |               ..     ..
%                          nnumnp |               xnumnp ynumnp];
%
%   ix:              = table of mesh connectivity information, specifying
%                      how nodes are attached to elements and how materials
%                      are assigned to elements; entries in the first nen
%                      columns correspond to the rows of NodeTable
%                      representing the nodes attached to element e;
%                      entries in the last nen+1 column are rows from MateT
%                      signifying the material properties assigned to
%                      element e; format of the table is as follows:
%                          Elements  |         n1    n2    n3    n4   mat
%                          e1        |  ix = [e1n1  e1n2  e1n3  e1n4 e1mat
%                          e2        |        e2n1  e2n2  e2n3  e2n4 e2mat
%                          ...       |         ..    ..    ..    ..   ..
%                          enumel    |        values for element numel   ];
%
%   MateT:           = table of mesh material properties for each distinct
%                      set of material properties; these sets are
%                      referenced by element e by setting the value of
%                      ix(e,nen+1) to the row number of the desired
%                      material set; format of the table is as follows:
%                          Materials  |           E   v   t
%                          mat1       |  MateT = [E1  v1  t1
%                          mat2       |           E2  v2  t2
%                          ...        |           ..  ..  ..];
%
%   BCLIndex:        = list of the number of boundary conditions and loads
%                      applied to the mesh; first entry is the number of
%                      prescribed displacements at nodes; second entry is
%                      the number of nodal forces
%
%   NodeBC:          = table of prescribed nodal displacement boundary
%                      conditions; it contains lists of nodes, the
%                      direction of the displacement prescribed (x=1, y=2),
%                      and the value of the displacement (set 0 for fixed
%                      boundary); the length of the table must match the
%                      entry in BCLIndex(1), otherwise an error will result
%                      if too few conditions are given or extra BCs will be
%                      ignored in the model input module;  format of the 
%                      table is as follows:
%                          BCs  |            nodes direction value
%                          bc1  |  NodeBC = [bc1n   bc1dir   bc1u
%                          bc2  |            bc2n   bc2dir   bc2u
%                          ...  |             ..     ..       .. ];
%
%   NodeLoad:        = table of prescribed nodal forces; it contains lists 
%                      of nodes, the direction of the force prescribed 
%                      (x=1, y=2), and the value of the force; the length 
%                      of the table must match the entry in BCLIndex(2), 
%                      otherwise an error will result if too few conditions
%                      are given or extra loads will be ignored in the 
%                      model input module; format of the table is as
%                      follows:
%                          Loads  |              nodes direction value
%                          P1     |  NodeLoad = [ P1n    P1dir    P1P
%                          P2     |               P2n    P2dir    P2P
%                          ...    |               ..     ..       .. ];
%
% The following numbering convention is used for 4-node quadrilateral
% elements:
%
%           4 -------------- 3       2
%           |                |       | \
%           |                |       |  \
%           |                |       |   \
%           |                |       |    \
%           |                |       |     \
%           |                |       |      \
%           1 -------------- 2       3-------1
%

format compact
clc

iel = 1; 
ndf = 3;
ndm = 3;

strainplot   = zeros(6,NbLoadSteps+1);
stressplot  = zeros(6,NbLoadSteps+1);
dispplot     = zeros(3,NbLoadSteps+1);
%qplot         = zeros(6,NbLoadSteps+1);
C1111_plot = zeros(1,NbLoadSteps+1);
C2222_plot = zeros(1,NbLoadSteps+1);
C1212_plot = zeros(1,NbLoadSteps+1);
res_count = 1;
fprintf(' Loadstep | NR-Iteration |     Abs. norm     |     Rel. Norm\n');

for step=1:NbLoadSteps
% for step=1:43    
%for alfa=1:12
%alfa=1;    
    if BCLIndex(1)>0
        NodeBC(1:BCLIndex(1),1:2)= NodeBC1(1:BCLIndex(1),1:2);
        NodeBC(1:BCLIndex(1),3)= NodeBC1(1:BCLIndex(1),3).*step./NbLoadSteps;
    end
    if BCLIndex(2)>0
        NodeLoad(1:BCLIndex(2),1:2)= NodeLoad1(1:BCLIndex(2),1:2);
        NodeLoad(1:BCLIndex(2),3)= NodeLoad1(1:BCLIndex(2),3).*step./NbLoadSteps;
    end
   
% Interpret Boundary Conditions and assign Loads; allocate dof's
assign_bc_load_data
Fext = Fd;
nneq = neq + nieq;

if step==1
    ModelDx         = zeros(neq, 1);
    %ep(1:6,1:8,1)  = zeros(6,8);
    %q(1:6,1:8,1)    = zeros(6,8);
    %alpha(1:8,1)    = zeros(8,1);
end

%Assemble Stiffness Routine
Fd = zeros(neq,1);
isw = 3;
FormFE
Fint=Fd;

%Solve Matrix System for FE Solution

SolveFE

R0=norm(Fdtilda);
Residual_log(res_count,1) = log10(1);
res_count=res_count+1;
fprintf('   %4d    |     %3d      |   %1.7e   |   %1.7e\n',...
            step,0,R0,1);
%R0=norm(Fext - Fint)
for iter=1:maxiterations
    Fd = zeros(neq,1);
    FormFE
    Fint=Fd;
    SolveFE
    R=norm(Fdtilda);
    Residual_log(res_count,1) = log10(R/R0);
    res_count=res_count+1;
    fprintf('   %4d    |     %3d      |   %1.7e   |   %1.7e\n',...
            step,iter,R,R/R0);
    if R<R0.*tol || R<10^-14 || iter == maxiterations
        iterations(step)=iter;
        break
    end
end
fprintf('\n');
Node_U_V = zeros(numnp,ndf);

for node = 1:numnp
    for dir = 1:ndf
        gDOF = NDOFT(node, dir);
        if gDOF <= neq
            Node_U_V(node, dir) = ModelDx(gDOF,1);
        else
            Node_U_V(node, dir) = ModelDc(gDOF - neq);
        end
    end
end

dispplot(1:3,step+1)   = Node_U_V(8, 1:3);

strainplot(1:6,step+1) = strain(1:6,3);
stressplot(1:6,step+1) = stress(1:6,3);
%qplot(1:6,step+1)       = q(1:6,2,step+1);
C1111_plot(1,step+1)    = C_ep(1,1,2,step+1);
C2222_plot(1,step+1)    = C_ep(2,2,2,step+1);
C1212_plot(1,step+1)    = C_ep(4,4,2,step+1);

if step == 1
    C1111_plot(1,1)    = C_ep(1,1,2,step+1);
    C2222_plot(1,1)    = C_ep(2,2,2,step+1);
    C1212_plot(1,1)    = C_ep(4,4,2,step+1);
end

if max(abs(strain),[],'all') >= limitingstrain
    break
end

end

figure (1)
plot(t,stressplot(1,:),'b','linewidth',3)
hold on
plot(t,stressplot(2,:),':r','linewidth',3)
plot(t,stressplot(4,:),'--g','linewidth',2)
plot(t,stressplot(3,:),'k','linewidth',2)
hold off
grid on
xlabel ('time (s)')
xlim ([0 1.2])
legend('\sigma_{11}','\sigma_{22}','\sigma_{12}','\sigma_{33}', 'location', 'northwest')

figure (2)
plot(strainplot(2,:),stressplot(1,:),'b','linewidth',3)
hold on
plot(strainplot(2,:),stressplot(2,:),':r','linewidth',3)
plot(strainplot(2,:),stressplot(4,:),'--g','linewidth',2)
plot(strainplot(2,:),stressplot(3,:),'k','linewidth',2)
hold off
grid on
xlabel ('E_{22}')
xlim ([0 0.25])
legend('\sigma_{11}','\sigma_{22}','\sigma_{12}','\sigma_{33}', 'location', 'northwest')

figure (3)
plot(strainplot(4,:),stressplot(1,:),'b','linewidth',3)
hold on
plot(strainplot(4,:),stressplot(2,:),':r','linewidth',3)
plot(strainplot(4,:),stressplot(4,:),'--g','linewidth',2)
plot(strainplot(4,:),stressplot(3,:),'k','linewidth',2)
hold off
grid on
xlabel ('E_{12}')
xlim ([0 0.25])
legend('\sigma_{11}','\sigma_{22}','\sigma_{12}','\sigma_{33}', 'location', 'northwest')


figure(4)
plot(Residual_log,'-o')
xlabel('iterations')
ylabel('log_{10}(Relative Residual)')


% Plot deformed configuration
factor=1;
NodeTable2 = NodeTable;
for i = 1:ndm
    NodeTable2(:,i) = NodeTable2(:,i) + Node_U_V(:,i)*factor;
end
plotModel(NodeTable2, ix, numel, nen, 5, 1, 1, 'Deformed Configuration (2D)', 'y', 'y')
plotModel(NodeTable, ix, numel, nen, 5, 1, 1, 'Deformed Configuration (2D)', 'n', 'n')