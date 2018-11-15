function [Match_Indexes,Best_Match_Indexes,T] = RANSAC_3D_Translation(SourceLocations,TargetLocations)
% ============================== About ====================================
% -------------------------------------------------------------------------
% Purpose: 3D Ransac
% Input: 
%   SourceLocations    : N*3,  The locations of the source features
%   TargetLocations    : N*3,  The locations of the target features
%
%  Output:
%   H                  : 3*3, Transformation Matrix
%   Match_Indexes      : N*2, Matched indexes   
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants  
% Northeastern University, USA
% =========================================================================
% -------------------------------------------------------------------------
MaxErrorDistance = 1;

Match_Indexes = [];
N=zeros(size(SourceLocations,2),1);
Translation=zeros(size(SourceLocations));
for i=1:size(SourceLocations,2)
    RandTargetLocations = TargetLocations(:,i);
    RandSourceLocations = SourceLocations(:,i);
    
    Translation(:,i) = RandTargetLocations - RandSourceLocations;
    TransformedPoints = TargetLocations - Translation(:,i);
    
    AllDistances = (sum((SourceLocations-TransformedPoints).^2,1)).^0.5;
    
   % CorrectNumbers = (AllDistances <= MaxErrorDistance);
    CorrectNumbers = find(AllDistances < MaxErrorDistance);
    N(i)=nnz(CorrectNumbers);
    
     if length(CorrectNumbers) > length(Match_Indexes)
         Match_Indexes = CorrectNumbers;
         T = Translation(:,i);
     end
end



% M=[2^0.5,3^0.5,5^0.5]*(Translation(:,N==max(N)));
% Match_Indexes=find(M==mode(M),1,'first');
Best_Match_Indexes=find(N==max(N),1,'first');
end

