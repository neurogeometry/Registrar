function NeighInd = findStackNeighbors(InputStackNumber,StackSizes_mm,StackPositions)
% ============================== About ====================================
% -------------------------------------------------------------------------
% Purpose: Find the nearest stacks of each stack
% Input: 
%   InputStackNumber    : 1*1,  The index of the stack to find neighbor
%                         stacks
%   StackSizes_mm       : N*3, (x,y,z)where N is the number of all stacks
%   StackPositions      : N*3, (x,y,z,)where N is the number of all stacks
%
%  Output:
%   NeighInd            : 1*N vector of all neighbors, 
%                         where N is the number of neighbors  
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants  
% Northeastern University, USA
% =========================================================================
% -------------------------------------------------------------------------
N=size(StackPositions,1);
StackCenters=StackPositions+StackSizes_mm./2;
StackCenters0=StackCenters(InputStackNumber,:);
StackSizes_mm0=StackSizes_mm(InputStackNumber,:);
ind=(abs(StackCenters-ones(N,1)*StackCenters0)<(StackSizes_mm+ones(N,1)*StackSizes_mm0)./2);
ind(InputStackNumber,1)=0;
NeighInd=find(sum(ind,2)==3);
end