function plotModel(NodeTable, ix, numPatches, nen, ModelID, subp, numsubp, plotname, interior, clfflag)
% plotModel(NodeTable, PatchTable, PatchIndex, PatchConn, U, V, W,
% numPatches, ModelID, geomref)
%
% Copyright (C) Arif Masud and Tim Truster
% UIUC
% 01/02/2009
%
%
%close all

%Set Face Colors
facecolors = [1 1 0     %yellow  U1 1
              0 0.8 0   %green   U2 2
              1 0 1     %magenta V1 3
              0.8 0 0   %red     V2 4
              0 1 1     %cyan    W1 5
              0 0 1];   %blue    W2 6

%Set Plot Resolution (number of "equi-distant" points on face);
numu = 10;
numv = 10;

%Select and clear figure
%Select and clear figure
figure(ModelID)
if clfflag == 'y'
clf(ModelID)
end
subplot(numsubp,1,subp)

if interior == 'n'
    linetyp = 'k:';
else
    linetyp = 'k-';
end

%Plot patch faces from ModelID


    for elem = 1:numPatches
    %     pause
    
        %Determine patch size parameters
        if nen < 5
            if nen == 3
                nel = 3;
                nelB = 3;
            elseif ix(elem,4) == 0
                nel = 3;
                nelB = 3;
            else
                nel = 4;
                nelB = 4;
            end
        else
            if nen == 6
                nel = 6;
                nelB = 3;
            elseif ix(elem,7) == 0
                nel = 6;
                nelB = 3;
            else
                nel = 9;
                nelB = 4;
            end
        end
        
        

            %Plot W Face
            PatchPoints = zeros(nel,2);
            
            %Extract face points
            for j = 1:nel
                    node = ix(elem,j);
                    for l = 1:2
                        PatchPoints(j,l) = NodeTable(node,l);
                    end %l
            end %j
          
            if nelB == 3
                len = 5;
                knots = [0 0.25 0.50 0.75 1.0 0 3/4/4 6/4/4 9/4/4 12/4/4 0 1/2/4 2/2/4 3/2/4 1/2 0 1/4/4 2/4/4 3/4/4 1/4 0 0 0 0 0
                         0 0 0 0 0 0.25 0.25 0.25 0.25 0.25 0.5  0.5  0.5  0.5  0.5 0.75 0.75 0.75 0.75 0.75 1.0  1.0  1.0  1.0  1.0]';
    %             plot2DLagransurfcont(PatchPoints, ContourValues, nel, knots, len, 'none', elem)
%                 plot2DLagransurfcont(PatchPoints, ContourValues, nel,
%                 knots, len, 'none', 0)
                if interior == 'y'
                plot2DLagransurf(PatchPoints, nel, 0, knots, len, 'none', elem)
                end
                hold on
                plot([PatchPoints(1,1) PatchPoints(2,1)],[PatchPoints(1,2) PatchPoints(2,2)],linetyp)
                plot([PatchPoints(2,1) PatchPoints(3,1)],[PatchPoints(2,2) PatchPoints(3,2)],linetyp)
                plot([PatchPoints(3,1) PatchPoints(1,1)],[PatchPoints(3,2) PatchPoints(1,2)],linetyp)
                hold off
            else
                len = 5;
                knots = [-1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0 -1.0 -0.5  0.0  0.5  1.0
                         -1.0 -1.0 -1.0 -1.0 -1.0 -0.5 -0.5 -0.5 -0.5 -0.5 -0.0 -0.0 -0.0 -0.0 -0.0  0.5  0.5  0.5  0.5  0.5  1.0  1.0  1.0  1.0  1.0]';
    %             if wf == 1 %NURBS
    %                 plot2DNURBSsurfcont(PatchPoints, ContourValues, UP, VP, nu, nv, pu, pv, nknotsU, nknotsV, Ulist, Vlist, 'none', n)
    %             else
    %                     plot2DLagransurfcont(PatchPoints, ContourValues, unum, vnum, numelU, numelV, pu, pv, nknotsU, nknotsV, 'none', ID)
                if interior == 'y'
                plot2DLagransurf(PatchPoints, nel, 0, knots, len, 'none', elem)
                end
                hold on
                plot([PatchPoints(1,1) PatchPoints(2,1)],[PatchPoints(1,2) PatchPoints(2,2)],linetyp)
                plot([PatchPoints(2,1) PatchPoints(3,1)],[PatchPoints(2,2) PatchPoints(3,2)],linetyp)
                plot([PatchPoints(3,1) PatchPoints(4,1)],[PatchPoints(3,2) PatchPoints(4,2)],linetyp)
                plot([PatchPoints(4,1) PatchPoints(1,1)],[PatchPoints(4,2) PatchPoints(1,2)],linetyp)
                hold off
            end

    end
title(plotname,'FontSize',14)
xlabel x
ylabel y
zlabel z