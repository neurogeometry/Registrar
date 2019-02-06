function [SubFeatures,SubFeaturesVectors] = OverlapRegion1(SourceID,TargetID,Source_seed_r_seed,Source_seed_FeatureVector,StackPositions_pixels,StackSizes_pixels,Pad)

SourceStackCenter=StackPositions_pixels(SourceID,:)+StackSizes_pixels(SourceID,:)./2;
TargetStackCenter=StackPositions_pixels(TargetID,:)+StackSizes_pixels(TargetID,:)./2;


PixelSizes = StackSizes_pixels(SourceID,:)./StackSizes_pixels(SourceID,1:3);

TargetDimX = StackSizes_pixels(TargetID,1);
TargetDimY = StackSizes_pixels(TargetID,2);
TargetDimZ = StackSizes_pixels(TargetID,3);

if SourceStackCenter(1) < TargetStackCenter(1)
    x1 = round((TargetStackCenter(1) - TargetDimX/2-StackPositions_pixels(SourceID,1))/PixelSizes(1))-Pad;
    x2 = round(StackSizes_pixels(SourceID,1)/PixelSizes(1));
elseif SourceStackCenter(1) > TargetStackCenter(1)
    x1 = 1;
    x2 =round((TargetStackCenter(1)+ TargetDimX/2 -StackPositions_pixels(SourceID,1))/PixelSizes(1))+Pad;
else
    x1 = 1;
    x2 = StackSizes_pixels(SourceID,1);
end


if SourceStackCenter(2) < TargetStackCenter(2)
    y1 = round((TargetStackCenter(2) - TargetDimY/2-StackPositions_pixels(SourceID,2))/PixelSizes(2))-Pad;
    y2 = round(StackSizes_pixels(SourceID,2)/PixelSizes(2));
elseif SourceStackCenter(2) > TargetStackCenter(2)
    y1 = 1;
    y2 =round((TargetStackCenter(2)+ TargetDimY/2 -StackPositions_pixels(SourceID,2))/PixelSizes(2))+Pad;
else
    y1 = 1;
    y2 = StackSizes_pixels(SourceID,2);
end



if SourceStackCenter(3) < TargetStackCenter(3)
    z1 = round((TargetStackCenter(3) - TargetDimZ/2-StackPositions_pixels(SourceID,3))/PixelSizes(3))-Pad;
    z2 = round(StackSizes_pixels(SourceID,3)/PixelSizes(3));
elseif SourceStackCenter(3) > TargetStackCenter(3)
    z1 = 1;
    z2 =round((TargetStackCenter(3)+ TargetDimZ/2 -StackPositions_pixels(SourceID,3))/PixelSizes(3))+Pad;
else
    z1 = 1;
    z2 = StackSizes_pixels(SourceID,3);
end

if StackSizes_pixels(1,3) >1
    ind = (Source_seed_r_seed(:,1)>x1 & Source_seed_r_seed(:,1)<x2 & Source_seed_r_seed(:,2)>y1 & Source_seed_r_seed(:,2)<y2 & Source_seed_r_seed(:,3)>z1 & Source_seed_r_seed(:,3)<z2);
else
    ind = (Source_seed_r_seed(:,1)>x1 & Source_seed_r_seed(:,1)<x2 & Source_seed_r_seed(:,2)>y1 & Source_seed_r_seed(:,2)<y2);
end
    
SubFeatures = Source_seed_r_seed(ind(:,1),:);
SubFeaturesVectors = Source_seed_FeatureVector(ind(:,1),:);

end