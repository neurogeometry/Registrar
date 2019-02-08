function All_overlaps=FindOverlaps(handles,StackPositions_pixels,StackSizes_pixels,StackList)
% ============================== About ====================================
% -------------------------------------------------------------------------
% Purpose: Find the nearest stacks of all stack
% Input:
%   StackPositions      : N*3, (x,y,z), where N is the number of all stacks
%   StackSizes_mm       : N*3, (x,y,z)where N is the number of all stacks
%
%  Output:
%   All_overlaps        : N*N matrix of overlapes
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
% Northeastern University, USA
% =========================================================================
% -------------------------------------------------------------------------
N=size(StackPositions_pixels,1);
StackCenters=StackPositions_pixels+StackSizes_pixels(:,1:3)./2;
count=0;
ii=zeros(13*N,1);
jj=zeros(13*N,1);
vv=zeros(13*N,1);
for InputStackNumber=1:N
    StackCenters0=StackCenters(InputStackNumber,:);
    StackSizes_mm0=StackSizes_pixels(InputStackNumber,:);
    NeighInd = findStackNeighbors(InputStackNumber,StackSizes_pixels,StackPositions_pixels(:,1:3));
    for j=1:length(NeighInd)
        if NeighInd(j)>InputStackNumber
            count=count+1;
            StackCenters1=StackCenters(NeighInd(j),:);
            StackSizes_mm1=StackSizes_pixels(NeighInd(j),:);
            ii(count)=InputStackNumber;
            jj(count)=NeighInd(j);
            vv(count) = prod(((StackSizes_mm1+StackSizes_mm0)./2)-abs(StackCenters0-StackCenters1));
        end
    end
end
All_overlaps = sparse(ii(1:count),jj(1:count),vv(1:count),N,N);

if handles.checkbox16.Value
    try
        DatasetMap();
        DatasetMapHandle=findobj(0,'Name','Stack Map');
        if ~isempty(DatasetMapHandle)
            cla(DatasetMapHandle,'reset')
            axis equal;
            axis off;
            view(3);
            DrawStackMap(StackPositions_pixels, StackSizes_pixels)
            for j=1:size(StackPositions_pixels,1)
                try
                    boxsize = [StackSizes_pixels(j,1),StackSizes_pixels(j,2),StackSizes_pixels(j,3)];
                    if ~isempty(StackList)
                        [PathStr,StackName]=fileparts(StackList{j,1});
                    else
                        StackName = '';
                    end
                    text(StackPositions_pixels(j,1)+boxsize(1)/2,StackPositions_pixels(j,2)+boxsize(2)/2,StackPositions_pixels(j,3)+boxsize(3)/2,[num2str(j),' ',StackName],'FontSize',12),hold on
                    drawnow();
                catch
                end
            end
        end
    catch
    end
end
