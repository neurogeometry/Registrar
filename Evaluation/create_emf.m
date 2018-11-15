%This script demonstrates some tricks to create an undistorted png image of
%the same size from imagesc.

function [] =create_emf(IM,AM,r)

outputsize_x = 1000;
outputsize_y = 1000;
figure,imshow(squeeze(max(IM(:,:,:),[],3)),[])
hold on
PlotAM_fast(AM,r)
h=findobj(gcf,'Type','Line');
% set(h,'Color',[1,0,0])
set(h,'LineWidth',1)
set(gcf,'Position',[0 0 outputsize_x outputsize_y],'Color',[0 0 0])
set(gca,'Position',[0 0 1 1])
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [outputsize_x outputsize_y]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 outputsize_x outputsize_y]);
print(gcf, '-dpng', 'xy-Projection.png');

