function plot2DLagransurfcont(P, Cont, nel, edge, knots, numu, gridlin, ID)
%(P, Cont, unum, vnum, numelU, numelV, pU, pV, numu, numv, Ulist, Vlist, gridlin, ID)
% plot3DNURBSsurf(Pw, U, V, nU, nV, pU, pV, numu, numv, facecol, gridlin)
%
% Copyright (C) Arif Masud and Tim Truster
% UIUC
% 06/22/2009
%
%

% SPL2D = zeros(numu,numv,2);
% CPL2D = zeros(numu,numv);
% 
% %Generate list of curve points
% 
% [SurPoint,CPL2D(1,1)] = LagrSurfacePointCont(nU,pU,U,nV,pV,V,P,Cont,-1,-1,2);
% for dir = 1:2
%     SPL2D(1,1,dir) = SurPoint(1,1,dir);
% end
% for i = 1:numelU
%     u = 2/unum*(1:unum) - ones(1,unum);
%     [SurPoint,CPL2D(i*unum-unum+2:i*unum+1,1)] = LagrSurfacePointCont(unum,pU,1,pV,P(i*pU-pU+1:i*pU+1,1,:),Cont(i*pU-pU+1:(i-1)*pU+1,1),u,-1,2);
%     for dir = 1:2
%         SPL2D(i*unum-unum+2:i*unum+1,1,dir) = SurPoint(:,:,dir);
%     end
%     for j = 1:numelV
%         v = 2/vnum*(1:vnum) - ones(1,vnum);
%         [SurPoint,CPL2D(i*unum-unum+2:i*unum+1,j*vnum-vnum+2,j*vnum+1)] = LagrSurfacePointCont(unum,pU,vnum,pV,P(i*pU-pU+1:i*pU+1,j*pV-pV+1,j*pV+1,:),Cont(i*pU-pU+1:i*pU+1,j*pV-pV+1,j*pV+1),u,v,2);
%         for dir = 1:2
%             SPL2D(i*unum-unum+2:i*unum+1,j*vnum-vnum+2:j*vnum+1,dir) = SurPoint(:,:,dir);
%         end
%     end
% end
% 
% %Plot the list of curve points, control points, and knots
% hold on
% surf(SPL2D(:,:,1),SPL2D(:,:,2),CPL2D,CPL2D,'EdgeColor',gridlin)
% uinc = round((numu-1)/2)+1;
% vinc = round((numv-1)/2)+1;
% text(SPL2D(uinc,vinc,1),SPL2D(uinc,vinc,2),CPL2D(uinc,vinc),num2str(ID))
% hold off

SPL2D = zeros(numu,numu,2);
CPL2D = zeros(numu,numu);

%Generate list of curve points
ind = 0;
for j = 1:numu
    for i = 1:numu
        ind = ind + 1;
        r = knots(ind,1);
        s = knots(ind,2);
        [SurPoint,CPL2D(i,j)] = LagrSurfacePointCont(P,Cont,nel,edge,r,s,2);
        for dir = 1:2
            SPL2D(i,j,dir) = SurPoint(dir);
        end
    end
end

%Plot the list of curve points, control points, and knots
hold on
surf(SPL2D(:,:,1),SPL2D(:,:,2),zeros(numu,numu),CPL2D,'EdgeColor',gridlin)
uinc = round((numu-1)/2)+1;
vinc = round((numu-1)/2)+1;
text(SPL2D(uinc,vinc,1),SPL2D(uinc,vinc,2),0,num2str(ID))
hold off