function [Match_Indexes] = RANSAC(SourceLocations,TargetLocations,TransformationValue,mu,NumRandPoints,MaxNumSamples,MaxNumMatches,MaxErrorDistance,InlierRatio,nxyz,affine)
NumAllMatches = size(TargetLocations,2);
if nchoosek(NumAllMatches,NumRandPoints)<=MaxNumSamples
    AllSamples=nchoosek(1:NumAllMatches,NumRandPoints);
else
    AllSamples=randi(NumAllMatches,[MaxNumSamples,NumRandPoints]);
end

NumCorrectHungarian = ceil(InlierRatio*NumAllMatches);
MaxNumMatches=max(MaxNumMatches,NumCorrectHungarian);
MaxNumMatches=min(MaxNumMatches,NumAllMatches);

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
        [~,L,b]=Optimal_Affine_Transform(RandSourceLocations,RandTargetLocations,mu);
        SourceLocations_affine=L*SourceLocations+b*ones(1,size(SourceLocations,2));
        AllDistances2 = sum((SourceLocations_affine-TargetLocations).^2,1);
    elseif TransformationValue == 4 % Non-Rigid

        [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(RandSourceLocations,RandTargetLocations,nxyz,affine,mu);
        [SourceLocations_nonRigid,~]=Perform_Bspline_Transform(SourceLocations,[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
        AllDistances2 = sum((SourceLocations_nonRigid-TargetLocations).^2,1);
    end
    
    CorrectNumbers = (AllDistances2 < MaxErrorDistance^2);
    i = i+1;
    if sum(CorrectNumbers) > length(Match_Indexes)
        Match_Indexes = find(CorrectNumbers);
    end
end
end