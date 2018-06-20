clear all
clc
addpath ('../functions');
Source_Stack_File = 'E:\SubTiling\00760\00760-ngc.0.tif';
IM_Source=ImportStack(char(Source_Stack_File)); 


points = [];
tic
for i=1:size(IM_Source,3)
tmp = detectMinEigenFeatures(IM_Source(:,:,i));
points = [points;[tmp.Location,ones(size(tmp.Location,1),1)*i]];
end
toc
imshow(max(IM_Source,[],3)); hold on;
plot(points(:,1),points(:,2),'g*');


points = [];
tic
for i=1:size(IM_Source,3)
tmp = detectMSERFeatures(IM_Source(:,:,i));
points = [points;[tmp.Location,ones(size(tmp.Location,1),1)*i]];
end
toc
imshow(max(IM_Source,[],3)); hold on;
plot(points(:,1),points(:,2),'g*');

I = max(IM_Source,[],3);
points = detectSURFFeatures(I,'MetricThreshold',50);
imshow(I); hold on;
plot(points.Location(:,1),points.Location(:,2),'g*');

points = detectKAZEFeatures(I);
imshow(I); hold on;
plot(points.Location(:,1),points.Location(:,2),'g*');

points = detectMSERFeatures(I);
imshow(I); hold on;
plot(points.Location(:,1),points.Location(:,2),'g*');





% for i=1:size(IM_Source,3)
% imwrite(IM_Source(:,:,i),['E:\JPEG2000Tiles\00760\00760-ngc_',num2str(i),'.jp2'],'jp2','Mode','lossy','CompressionRatio',20);
% end
% NEWStack = zeros(size(IM_Source));
% for i=1:size(IM_Source,3)
% NEWStack(:,:,i) = imread(['E:\JPEG2000Tiles\00760\00760-ngc_',num2str(i),'.jp2']);
% % figure(1);imshow(I);
% end
% 
% figure;imshow(max(NEWStack,[],3),[0 max(NEWStack(:))]);