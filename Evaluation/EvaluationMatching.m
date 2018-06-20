
    temp = sourceID;
        sourceID = targetID;
        targetID = temp;
    featuresFolder = '../data/features/';
    Source_seed = load([featuresFolder,char(StackList(sourceID,1)) '_seeds.mat' ]);
    Target_seed = load([featuresFolder,char(StackList(targetID,1)) '_seeds.mat' ]);
    [SourceSubFeatures,SourceSubFeaturesVectors] = OverlapRegion(sourceID,targetID,Source_seed);
    [TargetSubFeatures,TargetSubFeaturesVectors] = OverlapRegion(targetID,sourceID,Target_seed);
    
    tifFile_1 = dir([char(StackList(sourceID,2)) '/*.tif' ]);
    Source_Stack_Folder = [tifFile_1(1).folder,'/'];
    Source_Stack_File = tifFile_1(1).name;
    disp(['Reading ',Source_Stack_File]);
    IM_Source=ImportStack([Source_Stack_Folder,Source_Stack_File]);
    
    tifFile_2 = dir([char(StackList(targetID,2)) '/*.tif' ]);
    Target_Stack_Folder = [tifFile_2(1).folder,'/'];
    Target_Stack_File = tifFile_2(1).name;
    disp(['Reading ',Target_Stack_File]);
    IM_Target=ImportStack([Target_Stack_Folder,Target_Stack_File]);
    PixelSizes = StackSizes_mm(sourceID,:)./StackSizes_pixels(sourceID,1:3);
    Source_StackPositions = StackPositions(sourceID,:)./PixelSizes;
    Target_StackPositions = StackPositions(targetID,:)./PixelSizes;
    
    
    
    TempSource_StackPositions = Source_StackPositions([2,1,3]);
    TempTarget_StackPositions = Target_StackPositions([2,1,3]);
    Displacement = TempSource_StackPositions-TempTarget_StackPositions;
    if abs(Displacement(1)) > 0
        Direction='vertical';
        overlap = 0;
    elseif abs(Displacement(2)) > 0
        Direction='horizontal';
        overlap = 0;
    else
        Direction='vertical';
        overlap = 1;
    end
    Y1 = 1024;
    IM_source_max=max(IM_Source,[],3);
    IM_target_max=max(IM_Target,[],3);
    M_source=max(IM_Source(:));
    stitched = appendimages(IM_target_max,IM_source_max,Direction);
    figure('Position', [100 100 size(stitched,2) size(stitched,1)]);
    colormap('gray');
    figure(1);imshow(stitched,[0 M_source]);
    hold on
    loc1 = SourceSubFeatures;
    loc2 = TargetSubFeatures;
    
    loc1T = loc2 +[9,Y1-P1_P2_DX,0];
    plot(loc1T(:,2),loc1T(:,1),'r*');
    
    %     plot(loc2(:,2),loc2(:,1),'g*');
    loc1(:,2) = loc1(:,2)+Y1;
    plot(loc1(:,2),loc1(:,1),'c*');
    
    
    %    D = pdist2(loc1T,loc1,'seuclidean')
    load('species.mat');
    x= loc1T;
    y = loc1;
    
    
    
%     compute_nearest_neighbour(x,y)
    
    n_points = 68;%size(x,1)
%     figure
%  source = zeros(95,2);  
target = x;
source = y;

 K = 5;
%     source(:,:) = x(:,:);
   
target=target(:,[1,2]);
source=source(:,[1,2]);
[target_indices,dists] = knnsearch(target,source,'k',1,'distance','Minkowski','p',3)
    
 %   target_indices = compute_nearest_neighbour(source,target); 
    colors = hsv(n_points); 
  figure(1)
  hold on 
  for i = 1 :size(target_indices) 
      if target_indices(i)>0
     plot(source(i,2), source(i,1), 'o', 'Color', colors(i,:)) 
     
     plot(target(target_indices(i),2), target(target_indices(i),1), 's', ... 
        'Color', colors(i,:)) 
     plot( [ source(i,2) target(target_indices(i),2) ] , ... 
        [ source(i,1) target(target_indices(i),1) ], ... 
         'Color', colors(i,:)) 
      end
  end 
    % [indices,dists] = findNearestNeighbors(x,y,5);
%     P= x(:,:,:);
%     X=zeros(size(P));
%     X(1:68,:,:)= y;
%     I = nearestneighbour(X, P)
%     n = length(I);
%     plot(P(1,:), P(2, :), 'r.', X(1,:), X(2,:), 'b.', 'MarkerSize', 15)
%     hold on
%     
%     quiver(repmat(P(1), 1, n), repmat(P(2), 1, n), X(1, I) - P(1), X(2, I) - P(2), 0, 'k')
%     hold off
    



%    [n,d]=knnsearch(x,y,'k',5,'distance','minkowski','p',5);
% [ncb,dcb] = knnsearch(x,y,'k',10,...
%    'distance','chebychev');
% 
% 
%     gscatter(x(:,1),x(:,2),species)
% line(y(:,1),y(:,2),'marker','x','color','k',...
%    'markersize',10,'linewidth',2,'linestyle','none')
% line(x(n,1),x(n,2),'color',[.5 .5 .5],'marker','o',...
%    'linestyle','none','markersize',10)
% line(x(ncb,1),x(ncb,2),'color',[.5 .5 .5],'marker','p',...
%    'linestyle','none','markersize',10)
% legend('setosa','versicolor','virginica','query point',...
% 'minkowski','chebychev','Location','best')
    
