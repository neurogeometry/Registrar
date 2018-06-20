function [Match_Indexes] = RANSAC(SourceLocations,TargetLocations,TransformationValue)
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
if TransformationValue == 1 % Translation
    NumRandPoints = 1;
    MaxNumSamples=10^6;
    MaxNumMatches=inf;
    MaxErrorDistance = 2;
elseif TransformationValue == 2 % Rigid
    NumRandPoints = 2;
    MaxNumSamples=10^6;
    MaxNumMatches=inf;
    MaxErrorDistance = 2;
elseif TransformationValue == 3 % Affine
    NumRandPoints = 4;
    MaxNumSamples=10^6;
    MaxNumMatches=15;
    MaxErrorDistance = 2;
elseif TransformationValue == 4 % Non-Rigid
    NumRandPoints = 4;
%     N_L = 3;
    MaxNumSamples=10^3;
    MaxNumMatches=10;
    MaxErrorDistance = 4;
    Minimum = min(SourceLocations,[],2);
    Maximum = max(SourceLocations,[],2);
    nxyz = [512,512,60];
    Nxyz = ceil((Maximum-Minimum)./nxyz');
end

InlierRatio = 0.20;


NumAllMatches = size(TargetLocations,2);
if nchoosek(NumAllMatches,NumRandPoints)<=MaxNumSamples
    AllSamples=nchoosek(1:NumAllMatches,NumRandPoints);
else
    AllSamples=randi(NumAllMatches,[MaxNumSamples,NumRandPoints]);
end

NumCorrectHungarian = ceil(InlierRatio*NumAllMatches);
MaxNumMatches=max(MaxNumMatches,NumCorrectHungarian);
MaxNumMatches=min(MaxNumMatches,NumAllMatches);



% FailLimit = 100*nchoosek(NumAllMatches,NumRandPoints)/nchoosek(NumCorrectHungarian,NumRandPoints);
% if FailLimit > 10^5
%     disp('The Probability of Success in Ransac is to Low');
%     FailLimit = 10^5;
% end

i = 1;
Match_Indexes = [];
while  i <= size(AllSamples,1) && length(Match_Indexes) <= MaxNumMatches
    RandomSamples = AllSamples(i,:);
    RandTargetLocations = TargetLocations(:,RandomSamples);
    RandSourceLocations = SourceLocations(:,RandomSamples);
    
    if TransformationValue == 1 % Translation
        b=sum(RandTargetLocations-RandSourceLocations,2)./NumRandPoints;
        SourceLocations_Translation=SourceLocations+b*ones(1,size(SourceLocations,2));
        AllDistances2 = sum((SourceLocations_Translation-TargetLocations).^2,1);
    elseif TransformationValue == 2 % Rigid
        [L,b]=Optimal_Rigid_Transform(RandSourceLocations,RandTargetLocations);
        SourceLocations_Rigid=L*SourceLocations+b*ones(1,size(SourceLocations,2));
        AllDistances2 = sum((SourceLocations_Rigid-TargetLocations).^2,1);
    elseif TransformationValue == 3 % Affine
        [L,b]=Optimal_Affine_Transform(RandSourceLocations,RandTargetLocations);
        SourceLocations_affine=L*SourceLocations+b*ones(1,size(SourceLocations,2));
        AllDistances2 = sum((SourceLocations_affine-TargetLocations).^2,1);
    elseif TransformationValue == 4 % Non-Rigid
%         Min=min([min(RandSourceLocations,[],2),min(TargetLocations,[],2)],[],2);
%         Max=max([max(RandSourceLocations,[],2),max(TargetLocations,[],2)],[],2);
%         [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(RandSourceLocations,RandTargetLocations,N_L,Min,Max);
%         SourceLocations_nonRigid=Perform_Nonrigid_Transform(SourceLocations,XYZlmn,N_L,Min,Max);
%         AllDistances2 = sum((SourceLocations_nonRigid-TargetLocations).^2,1);

        [~,Cxyz]=Optimal_Bspline_Transform(RandSourceLocations-Minimum*ones(1,NumRandPoints),RandTargetLocations-Minimum*ones(1,NumRandPoints),Nxyz,nxyz);
        [SourceLocations_nonRigid,~]=Perform_Bspline_Transform(SourceLocations-Minimum*ones(1,size(SourceLocations,2)),Cxyz,Nxyz,nxyz);
        AllDistances2 = sum((SourceLocations_nonRigid+Minimum*ones(1,size(SourceLocations,2))-TargetLocations).^2,1);
    end
    
    CorrectNumbers = (AllDistances2 < MaxErrorDistance^2);
    i = i+1;
    if sum(CorrectNumbers) > length(Match_Indexes)
        Match_Indexes = find(CorrectNumbers);
    end
end

end