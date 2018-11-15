function [H, Match_Indexes] = RANSAC_3D(SourceLocations,TargetLocations)
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
NumRandPoints = 4;
MaxErrorDistance = 1;
InlierRatio = 0.15;
CorrectHungarianRatio = 0.3;


NumAllMatches = size(TargetLocations,2);
NumCorrectNeeded = round(InlierRatio*NumAllMatches);
NumCorrectHungarian = round(CorrectHungarianRatio*NumAllMatches);
 
FailLimit = 100*nchoosek(NumAllMatches,NumRandPoints)/nchoosek(NumCorrectHungarian,NumRandPoints);
if FailLimit > 10^6
    disp('The Probability of Success in Ransac is to Low');
    FailLimit = 10^6;
end

CorrectNumbers = [];
i = 1;
Match_Indexes = [];
H = [];
while  length(CorrectNumbers) < NumCorrectNeeded && i < FailLimit
    RandomSamples = randi([1 NumAllMatches],1,NumRandPoints);
    if size(SourceLocations,1) == 2
        H = Transformation2D(TargetLocations(:,RandomSamples),SourceLocations(:,RandomSamples));
        AllDistances = calcDist2D(H,TargetLocations,SourceLocations);
    else  
        [H,AllDistances] = Transformation3D(TargetLocations,SourceLocations,RandomSamples);
        %[AllDistances] = ShiftTransformation3D(TargetLocations,SourceLocations,RandomSamples);
        %AllDistances = calcDist3D(H,TargetLocations,SourceLocations);
    end
    
    CorrectNumbers = find(AllDistances < MaxErrorDistance);
    i = i+1;
    if length(CorrectNumbers) > length(Match_Indexes)
        Match_Indexes = CorrectNumbers;
    end
       
end

end




