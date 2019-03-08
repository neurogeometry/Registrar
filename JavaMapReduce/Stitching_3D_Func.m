function [MatchLocations,MatchLocationsHang] = Stitching_3D_Func(SourceID,TargetID,StackPositions_pixels,StackSizes_pixels,Source_seed_r_seed,Source_seed_FeatureVector,Target_seed_r_seed,Target_seed_FeatureVector,TransformationValue)

Source_StackPositions = StackPositions_pixels(SourceID,:);
Target_StackPositions = StackPositions_pixels(TargetID,:);
[Source_seed,SourceFeatures] = OverlapRegion1(SourceID,TargetID,Source_seed_r_seed,Source_seed_FeatureVector,StackPositions_pixels,StackSizes_pixels,0);
[Target_seed,TargetFeatures] = OverlapRegion1(TargetID,SourceID,Target_seed_r_seed,Target_seed_FeatureVector,StackPositions_pixels,StackSizes_pixels,0);
        


paramsFMDT = 60;
paramsFM_C1 = [];
paramsFM_C2 = [];
mu = 1024;
Displacement = Source_StackPositions-Target_StackPositions;

hungInput = zeros(size(Source_seed,1),size(Target_seed,1));


    nTargetFeatures=bsxfun(@rdivide,TargetFeatures,mean(TargetFeatures,2));
    nSourceFeatures=bsxfun(@rdivide,SourceFeatures,mean(SourceFeatures,2));
    for i=1:size(SourceFeatures,2)
        hungInput=hungInput+abs(bsxfun(@minus,nTargetFeatures(:,i)',nSourceFeatures(:,i)));
    end
    hungInput=hungInput./size(SourceFeatures,2);
    
    TempDisp=zeros(size(SourceFeatures,1),size(TargetFeatures,1));
    for i=1:size(Target_seed,2)
        TempDisp=TempDisp+(bsxfun(@minus,Target_seed(:,i)',Source_seed(:,i))).^2;
    end
    Dthr=(sum(Displacement.^2))^0.5;
    hungInput(isnan(hungInput)) = 10^12;
    hungInput(abs(TempDisp.^0.5-Dthr)>paramsFMDT)=10^12;


Am = Hungarian_fast(hungInput,paramsFM_C1,paramsFM_C2);
Am(hungInput==10^12)=0;
[idx1,idx2]=find(Am);

x_Source = Source_seed(idx1,1);
x_Target = Target_seed(idx2,1);
y_Source = Source_seed(idx1,2);
y_Target = Target_seed(idx2,2);
z_Source = Source_seed(idx1,3);
z_Target = Target_seed(idx2,3);


matchLoc_Target = [x_Target,y_Target,z_Target];
matchLoc_Source = [x_Source,y_Source,z_Source];


Match_Indexes=[];
Global_Matched_Source = matchLoc_Source'+(ones(size(matchLoc_Source,1),1)*Source_StackPositions)'-1;
Global_Matched_Target = matchLoc_Target'+(ones(size(matchLoc_Source,1),1)*Target_StackPositions)'-1;

MatchLocationsHang = [Global_Matched_Source',Global_Matched_Target'];

try
    Match_Indexes = RANSAC(Global_Matched_Source,Global_Matched_Target,TransformationValue,mu);
catch
    
end

MatchLocations=NaN(length(Match_Indexes),6);
if ~isempty(Match_Indexes)
    for i = 1:length(Match_Indexes)
        MatchLocations(i,:) = [matchLoc_Source(Match_Indexes(i),:),matchLoc_Target(Match_Indexes(i),:)];
    end
end

end