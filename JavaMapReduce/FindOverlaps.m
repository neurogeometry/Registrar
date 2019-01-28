function All_overlaps=FindOverlaps(StackPositions_pixels,StackSizes_pixels)

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