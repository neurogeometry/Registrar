function [Match_Indexes] = RANSAC_Translation(SourceLocations,TargetLocations)
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
NumRandPoints = 1;
MaxErrorDistance = 2;
InlierRatio = 0.20;

NumAllMatches = size(TargetLocations,2);
NumCorrectHungarian = ceil(InlierRatio*NumAllMatches);
FailLimit = 10*nchoosek(NumAllMatches,NumRandPoints)/nchoosek(NumCorrectHungarian,NumRandPoints);%10000;%
if FailLimit < 1000
    FailLimit = 1000;
elseif FailLimit > 10^5
    disp('The Probability of Success in Ransac is to Low');
    FailLimit = 10^5;
end

i = 1;
Match_Indexes = [];
while  i < FailLimit
    RandomSamples = randi([1 NumAllMatches],1,NumRandPoints);
    
    RandTargetLocations = TargetLocations(:,RandomSamples);
    RandSourceLocations = SourceLocations(:,RandomSamples);
    
%     [L,b]=Optimal_Affine_Transform(RandTargetLocations,RandSourceLocations);
    b=Optimal_Translation_Transform(RandSourceLocations,RandTargetLocations);
    SourceLocations_Translation=SourceLocations+b*ones(1,size(SourceLocations,2));
    AllDistances = sum((SourceLocations_Translation-TargetLocations).^2,1).^0.5;
    
    CorrectNumbers = find(AllDistances < MaxErrorDistance);
    i = i+1;
    if length(CorrectNumbers) > length(Match_Indexes)
        Match_Indexes = CorrectNumbers;
    end
       
end
% i
end




