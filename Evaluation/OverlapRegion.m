function [SubFeatures,SubFeaturesVectors] = OverlapRegion(TargetID,SourceID,Source_seed)
load('../data/StackData.mat');
SourceStackCenter=StackPositions(SourceID,:)+StackSizes_mm(SourceID,:)./2;
TargetStackCenter=StackPositions(TargetID,:)+StackSizes_mm(TargetID,:)./2;
Pad = -40;

PixelSizes = StackSizes_mm(SourceID,:)./StackSizes_pixels(SourceID,1:3);
% SourceDimX = StackSizes_mm(SourceID,1);
% SourceDimY = StackSizes_mm(SourceID,2);
% SourceDimZ = StackSizes_mm(SourceID,3);

TargetDimX = StackSizes_mm(TargetID,1);
TargetDimY = StackSizes_mm(TargetID,2);
TargetDimZ = StackSizes_mm(TargetID,3);

if SourceStackCenter(1) < TargetStackCenter(1)
    x1 = round((TargetStackCenter(1) - TargetDimX/2-StackPositions(SourceID,1))/PixelSizes(1))-Pad;
    x2 = round(StackSizes_mm(SourceID,1)/PixelSizes(1));
elseif SourceStackCenter(1) > TargetStackCenter(1)
    x1 = 1;
    x2 =round((TargetStackCenter(1)+ TargetDimX/2 -StackPositions(SourceID,1))/PixelSizes(1))+Pad;
else
    x1 = 1;
    x2 = StackSizes_pixels(SourceID,1);
end


if SourceStackCenter(2) < TargetStackCenter(2)
    y1 = round((TargetStackCenter(2) - TargetDimY/2-StackPositions(SourceID,2))/PixelSizes(2))-Pad;
    y2 = round(StackSizes_mm(SourceID,2)/PixelSizes(2));
elseif SourceStackCenter(2) > TargetStackCenter(2)
    y1 = 1;
    y2 =round((TargetStackCenter(2)+ TargetDimY/2 -StackPositions(SourceID,2))/PixelSizes(2))+Pad;
else
    y1 = 1;
    y2 = StackSizes_pixels(SourceID,2);
end



if SourceStackCenter(3) < TargetStackCenter(3)
    z1 = round((TargetStackCenter(3) - TargetDimZ/2-StackPositions(SourceID,3))/PixelSizes(3))-Pad;
    z2 = round(StackSizes_mm(SourceID,3)/PixelSizes(3));
elseif SourceStackCenter(3) > TargetStackCenter(3)
    z1 = 1;
    z2 =round((TargetStackCenter(3)+ TargetDimZ/2 -StackPositions(SourceID,3))/PixelSizes(3))+Pad;
else
    z1 = 1;
    z2 = StackSizes_pixels(SourceID,3);
end

% StackSizes_pixels(SourceID,1:3);
% StackSizes_pixels(TargetID,1:3);

ind = (Source_seed.r_seed(:,1)>y1 & Source_seed.r_seed(:,1)<y2 & Source_seed.r_seed(:,2)>x1 & Source_seed.r_seed(:,2)<x2 & Source_seed.r_seed(:,3)>z1 & Source_seed.r_seed(:,3)<z2);


SubFeatures = Source_seed.r_seed(ind,:);
SubFeaturesVectors = Source_seed.FeatureVector(ind,:);






end