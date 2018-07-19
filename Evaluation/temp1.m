clear all;
close all;
clc;

sourceID = 1;
targetID = 2;
TraceNum = 2;
mu = 1040;
affine = 1;
ppm = 1;
N = 1; %Image Numbers
addpath('NeuronTracerV20');
csvPath = '..\..\RegistrationEvaluation\TimeLapse_Holtmaat_StackList.csv';
StackList = table2cell(readtable(csvPath,'Delimiter',','));
GTpath = '..\..\RegistrationEvaluation\TimeLaps\';
Functionspath = '..\';
resultpath = '..\..\RegistrationEvaluation\Results-TimeLapse_Holtmaat_StackList\';
addpath([Functionspath,'Functions']);
% Laod all matched Boutons
load([GTpath,'Matches\DL083-001-Matches.mat']);


T_Names = {'B','C','D','E','F','G','H','I','J','K','L','M','N'};
T_Days = 1:4:41;

% calculate Nonrigid Transform
%  FeaturePositions_NR = load([resultpath,'MatchedPoints_Non-Rigid_mu0.mat']);

%FeaturePositions_NR = load('MatchedPoints_Non-Rigid_mu1020.mat');
% Matched_hung = FeaturePositions_NR.Matched_hung;

FeaturePositions_NR = load('MatchedPoints_Non-Rigid_mu1020_new.mat');
Matched_hung = FeaturePositions_NR.Matched_hung;

Matched = FeaturePositions_NR.Matched;

% Matched_hung = FeaturePositions_NR.Matched_hung;
Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';
nxyz = [256;256;156];


%Use Boutons
Boutons = struct2cell(Matches);
B1 = Boutons{sourceID,1};
SourcePoints = B1.r1;
TargetPoints = B1.r2;
BoutonsBefore = mean(sum((TargetPoints-SourcePoints).^2,2).^0.5)
% k=0;
% for i = 1:12
% B1 = Boutons{i,1};
% k = k + size(B1.r1,1);
% end

% Get all Traces
fname_First = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
fname_First={fname_First.name}';
fname_Second = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
fname_Second={fname_Second.name}';
color = rand(size(fname_First,1),3);

% Laod images
Source_Stack_File = char(StackList(sourceID,1));
Target_Stack_File = char(StackList(targetID,1));
IM_Source=ImportStack(char(Source_Stack_File));
IM_Source = uint8(double(IM_Source)./double(max(IM_Source(:))).*255);
IM_Target=ImportStack(char(Target_Stack_File));
IM_Target = uint8(double(IM_Target)./double(max(IM_Target(:))).*255);
IM_source_max=max(IM_Source,[],3);
IM_target_max=max(IM_Target,[],3);

%% Show All Traces  ----------------- Fig 1. A
% N = N + 1;
% figure(N);imshow(IM_source_max,[0 50]);
% N = N + 1;
% figure(N);imshow(IM_target_max,[0 50]);
% for i=1:size(fname_First,1)
%     % Using Trace
%     
%     sourcePath = [GTpath,'Matches\Traces\',fname_First{i}];
%     targetPath = [GTpath,'Matches\Traces\',fname_Second{i}];
%     
%     [AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
%     [AM_Target,r_Target,R_Target]=swc2AM(targetPath);
%     [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
%     [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
%     
%     figure(N-1);hold on; PlotAM(AM_Source,r_Source,color(i,:))
%     text(r_Source(round(size(r_Source,1)/3),2),r_Source(round(size(r_Source,1)/3),1),['Axon #',num2str(i)],'Color',color(i,:))
%     figure(N);hold on; PlotAM(AM_Target,r_Target,color(i,:))
%     text(r_Target(round(size(r_Target,1)/3),2),r_Target(round(size(r_Target,1)/3),1),['Axon #',num2str(i)],'Color',color(i,:))
% 
% end

% %% Show One Trace  ----------------- Fig 1. B
% TraceNum = 2;
% % for TraceNum=1:22
% sourcePath = [GTpath,'Matches\Traces\',fname_First{TraceNum}];
% targetPath = [GTpath,'Matches\Traces\',fname_Second{TraceNum}];
% [AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
% [AM_Target,r_Target,R_Target]=swc2AM(targetPath);
% [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
% [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
% 
% [KT]=FastMarchingTube(size(IM_Source),r_Source,7,[1,1,1]);
% N = N+1;
% figure(N);imshow(max(uint8(KT).*IM_Source,[],3),[0 30])
% figure(N);hold on; PlotAM(AM_Source,r_Source,color(TraceNum,:))
% % figure(N);hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'or')
% % text(r_Source(round(size(r_Source,1)/3),2),r_Source(round(size(r_Source,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(0)],'Color',color(TraceNum,:))
% 
% [KT]=FastMarchingTube(size(IM_Target),r_Target,7,[1,1,1]);
% N = N+1;
% figure(N);imshow(max(uint8(KT).*IM_Target,[],3),[0 30])
% figure(N);hold on; PlotAM(AM_Target,r_Target,color(TraceNum,:))
% % figure(N);hold on; plot(TargetPoints(:,2),TargetPoints(:,1),'or')
% % text(r_Target(round(size(r_Target,1)/3),2),r_Target(round(size(r_Target,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(4)],'Color',color(TraceNum,:))
% % end
%% Show Overlaped Images Fig. 1. C
% N = N+1;
% figure(N),imshowpair(IM_source_max,IM_target_max,'Scaling','independent')

%% Show One Trace  ----------------- Fig 1. D
% 
% for i = 1: size(StackList,1)
%     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{i},'001-A0*']);
%     fname={fname.name}';
%     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%     [AM,r,R]=swc2AM(Path);
%     [AM,r,~] = AdjustPPM(AM,r,R,ppm);
%     if i ==1
%     Stack_File = char(StackList(i,1));
%     IM=ImportStack(char(Stack_File));
%     IM = uint8(double(IM)./double(max(IM(:))).*255);
%     IM_max=max(IM,[],3);
% 
%     [KT]=FastMarchingTube(size(IM),r,3,[1,1,1]);
% %     N = N+1;
%     figure(N);imshow(max(uint8(KT).*IM,[],3),[0 20])
%     end
%     figure(N); hold on; PlotAM(AM,r,color(TraceNum+i-1,:))
% %     text(r(round(size(r,1)/3),2),r(round(size(r,1)/3),1),['Day ',num2str(T_Days(i))],'Color',color(TraceNum+i-1,:))
% %     text(r(round(size(r,1)/3),2),r(round(size(r,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(T_Days(i))],'Color',color(TraceNum,:))
%     if i == size(StackList,1)
%         B1 = Boutons{i-1,1};
%         SourcePoints = B1.r2;
%         
%     else
%         B1 = Boutons{i,1};
%         SourcePoints = B1.r1;
%         
%     end
% %     figure(N);hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'or')
% end

% % Matching Result  ----------------- Fig 2.
% stitched = appendimages(IM_target_max,IM_source_max,'horizontal');
% %                   ----------------- Fig 2. A
% N = N + 1;
% matchesStep = 3;
% figure(N);imshow(stitched,[0 20])
% figure(N);hold on;plot(Matched_hung{sourceID,targetID}(1:matchesStep:end,5),Matched_hung{sourceID,targetID}(1:matchesStep:end,4),'*r');
% figure(N);hold on;plot(Matched_hung{sourceID,targetID}(1:matchesStep:end,2)+size(IM_target_max,1)+45,Matched_hung{sourceID,targetID}(1:matchesStep:end,1),'*r');
% figure(N);hold on;
% for i = 1:matchesStep:length(Matched_hung{sourceID,targetID})
% line([Matched_hung{sourceID,targetID}(i,5) Matched_hung{sourceID,targetID}(i,2)+size(IM_target_max,1)+45], ...
%                     [Matched_hung{sourceID,targetID}(i,4) Matched_hung{sourceID,targetID}(i,1)], 'Color', rand(1,3));
% end
% 
% %                   ----------------- Fig 2. B
% N = N + 1;
% figure(N);imshow(stitched,[0 20])
% figure(N);hold on;plot(Matched{sourceID,targetID}(:,5),Matched{sourceID,targetID}(:,4),'*r');
% figure(N);hold on;plot(Matched{sourceID,targetID}(:,2)+size(IM_target_max,1)+45,Matched{sourceID,targetID}(:,1),'*r');
% figure(N);hold on;
% for i = 1: length(Matched{sourceID,targetID})
% line([Matched{sourceID,targetID}(i,5) Matched{sourceID,targetID}(i,2)+size(IM_target_max,1)+45], ...
%                     [Matched{sourceID,targetID}(i,4) Matched{sourceID,targetID}(i,1)], 'Color', rand(1,3));
% end
% %                   ----------------- Fig 2. C
% [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
% nx=10; ny=10; nz=10;
% Min=min(Global_Matched_Source,[],2);
% Max=max(Global_Matched_Source,[],2);
% [xx,yy,zz]=ndgrid(1:nx,1:ny,1:nz);
% r_grid=[Min(1)+(Max(1)-Min(1)).*(xx(:)-1)./(nx-1),Min(2)+(Max(2)-Min(2)).*(yy(:)-1)./(ny-1),Min(3)+(Max(3)-Min(3)).*(zz(:)-1)./(nz-1)]';
% AM_grid=(abs(bsxfun(@minus,xx(:),xx(:)'))+abs(bsxfun(@minus,yy(:),yy(:)'))+abs(bsxfun(@minus,zz(:),zz(:)')))==1;
% N = N +1;figure(N);
% subplot(1,2,1)
% PlotAM(AM_grid,r_grid','k')
% title('Original')
% [r_grid_aligned,~,~]=Perform_Bspline_Transform(r_grid,[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
% subplot(1,2,2)
% PlotAM(AM_grid,r_grid_aligned','k')
% title('Aligned')

%% Validation Result  ----------------- Fig 3.A
% [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
% 
% [IM_Target_NR,StackPosition_prime,~]=Perform_Bspline_Transform(IM_Target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
% IM_Target_NR_max=max(IM_Target_NR,[],3);
% MIN=min([1;1],StackPosition_prime(1:2));
% MAX=max(size(IM_source_max)',size(IM_target_max)'+StackPosition_prime(1:2)-1);
% temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');
% 
% IM_Target_NR_max_P=temp;
% IM_Target_NR_max_P(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_Target_NR_max,1),...
%     StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_Target_NR_max,2))=IM_Target_NR_max;
% IM_source_max_P=temp;
% IM_source_max_P(2-MIN(1):1-MIN(1)+size(IM_source_max,1),...
%     2-MIN(2):1-MIN(2)+size(IM_source_max,2))=IM_source_max;
% N = N +1;
% figure(N),imshowpair(IM_Target_NR_max_P,IM_source_max_P,'Scaling','independent')

%                      ----------------- Fig 3.B
% L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
% for sourceID = 1: size(StackList,1)-1
%     
%     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
%     fname={fname.name}';
%     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%     [AM,r,R]=swc2AM(Path);
%     [AM,r,~] = AdjustPPM(AM,r,R,ppm);
%     
%      if sourceID ==1
%         Stack_File = char(StackList(sourceID,1));
%     IM=ImportStack(char(Stack_File));
%     IM = uint8(double(IM)./double(max(IM(:))).*255);
%     IM_max=max(IM,[],3);
%     [KT]=FastMarchingTube(size(IM),r,3,[1,1,1]);
%     N = N+1;
%     figure(N);imshow(max(uint8(KT).*IM,[],3),[0 20])
%     figure(N);hold on; PlotAM(AM,r,color(TraceNum,:))
% %     text(r(round(size(r,1)/3),2),r(round(size(r,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(T_Days(i))],'Color',color(TraceNum,:))
% %     if i == size(StackList,1)
% %         B1 = Boutons{i-1,1};
% %         SourcePoints = B1.r2;  
% %     else
% %         B1 = Boutons{i,1};
% %         SourcePoints = B1.r1;
% %     end
% %     figure(N);hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'or')
%      end
%    
%     targetID = sourceID + 1;
%     Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
%     Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';
%     [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
%     
%     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
%     fname={fname.name}';
%     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%     [AM,r,R]=swc2AM(Path);
%     [AM,r,~] = AdjustPPM(AM,r,R,ppm);
%     Points_NR_temp = r';
%     for j=size(b,2):-1:1
%     [Points_NR_temp,~]=Perform_Bspline_Transform(Points_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
%     end
%     hold on; PlotAM(AM,Points_NR_temp',color(TraceNum+sourceID,:))
% 
% %     [KT]=FastMarchingTube(size(IM),r,3,[1,1,1]);
% %     Target_Stack_File = char(StackList(targetID,1));
% %     IM_Target=ImportStack(char(Target_Stack_File));
% %     IM_Target = uint8(double(IM_Target)./double(max(IM_Target(:))).*255);
% %     IM_Target = uint8(KT).*IM_Target;
% %     [IM_Target_NR,StackPosition_prime,~]=Perform_Bspline_Transform(IM_Target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
% %     IM_Target_NR_max=max(IM_Target_NR,[],3);
% %     MIN=min([1;1],StackPosition_prime(1:2));
% %     MAX=max(size(IM_source_max)',size(IM_target_max)'+StackPosition_prime(1:2)-1);
% %     temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');
% %     IM_Target_NR_max_P=temp;
% %     IM_Target_NR_max_P(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_Target_NR_max,1),...
% %         StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_Target_NR_max,2))=IM_Target_NR_max;
% %     IM_source_max_P=temp;
% %     IM_source_max_P(2-MIN(1):1-MIN(1)+size(IM_source_max,1),...
% %         2-MIN(2):1-MIN(2)+size(IM_source_max,2))=IM_source_max;
% %     N = N +1;
% %     figure(N),imshowpair(IM_Target_NR_max_P,IM_source_max_P,'Scaling','independent')
% 
% end

%                      ----------------- Fig 3.C
%From Evaluate_TimeLapse.m file
%                      ----------------- Fig 3.D
%From Evaluate_TimeLapse.m file


%%                     ----------------- Fig 4.A
% IMAll = uint8(zeros(1024,1024));
% TR = 0;
% AddSpace = 40;
% for ID = 1: size(StackList,1)
%     sourceID = ID;
%     targetID = ID + 1;
%     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
%     fname={fname.name}';
%     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%     [AM_source,r_source,R]=swc2AM(Path);
%     [AM_source,r_source,~] = AdjustPPM(AM_source,r_source,R,ppm);
%     
%     
%     
%         Stack_File = char(StackList(sourceID,1));
%         IM=ImportStack(char(Stack_File));
%         IM = uint8(double(IM)./double(max(IM(:))).*255);
%         IM_max=max(IM,[],3);
%     
%         [KT]=FastMarchingTube(size(IM),r_source,3,[1,1,1]);
%         IMmax =  max(uint8(KT).*IM,[],3);
%         IMmax = imtranslate(IMmax,[0, TR]);
%         IMAll = max(IMAll,IMmax);
%     
%     if ID < size(StackList,1)
%         
%         fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
%         fname={fname.name}';
%         Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%         [AM_target,r_target,R]=swc2AM(Path);
%         [AM_target,r_target,~] = AdjustPPM(AM_target,r_target,R,ppm);
%         
%         B1 = Boutons{sourceID,1};
%         SourcePoints = B1.r1;
%         TargetPoints = B1.r2;
%         
%         
%         [k,f] = dsearchn(r_source,SourcePoints);
%         
%         s = 1;
%         SourcePoints_onTrace=[];
%         for j=1:size(k)
%             if f(j)< 3
%                 SourcePoints_onTrace(s,:) = SourcePoints(j,:);
%                 s = s +1;
%             end
%         end
%         figure(N); hold on;
%         %     if ID==1
%         %     figure(N); hold on; PlotAM(AM_source,r_source,'c')
%         %
%         %     end
%         
%         plot3(SourcePoints_onTrace(:,2),SourcePoints_onTrace(:,1)+TR,ones(1,length(SourcePoints_onTrace(:,2))),'or')
%         
%         TR = TR - AddSpace;
%         [k,f] = dsearchn(r_target,TargetPoints);
%         TargetPoints(:,1) = TargetPoints(:,1)+TR;
%         r_target(:,1)=r_target(:,1)+TR;
%         %     PlotAM(AM_target,r_target,'c')
%         s = 1;
%         TargetPoints_onTrace=[];
%         for j=1:size(k)
%             if f(j)< 3
%                 TargetPoints_onTrace(s,:) = TargetPoints(j,:);
%                 s = s +1;
%             end
%         end
%         plot3(TargetPoints_onTrace(:,2),TargetPoints_onTrace(:,1),ones(1,length(TargetPoints_onTrace(:,2))),'or')
%         
%         if ID==1
%             line([SourcePoints_onTrace(:,2) TargetPoints_onTrace(:,2)]', ...
%                 [SourcePoints_onTrace(:,1) TargetPoints_onTrace(:,1)]',[1;1]*ones(1,length(TargetPoints_onTrace(:,1))), 'Color', 'c');
%         else
%             line([SourcePoints_onTrace(:,2) TargetPoints_onTrace(:,2)]', ...
%                 [SourcePoints_onTrace(:,1)+TR+AddSpace TargetPoints_onTrace(:,1)]',[1;1]*ones(1,length(TargetPoints_onTrace(:,1))), 'Color', 'c');
%         end
%         drawnow
%     end
%     
% end
% figure(N);imshow(IMAll,[0 20])

%%                     ----------------- Fig 4.B

L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
IMAll = uint8(zeros(1024,1024));
TR = 0;
AddSpace = 40;
for ID = 1: size(StackList,1)
    sourceID = ID;
    targetID = ID + 1;
    
    Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';
    [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
    
    fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
    fname={fname.name}';
    Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
    [AM_source,r_source,R]=swc2AM(Path);
    [AM_source,r_source,~] = AdjustPPM(AM_source,r_source,R,ppm);
    
    Stack_File = char(StackList(sourceID,1));
    IM=ImportStack(char(Stack_File));
    IM = uint8(double(IM)./double(max(IM(:))).*255);
    IM_max=max(IM,[],3);
    
    [KT]=FastMarchingTube(size(IM),r_source,3,[1,1,1]);
    IM =  uint8(KT).*IM;
    IMind=find(KT);
    [IMx,IMy,IMz]=ind2sub(size(IM),IMind);
    IMxyz=[IMx,IMy,IMz]';
    Int=IM(IMind);
    
    for j=ID-1:-1:1
        StackPosition_prime=[1;1;1];
        [IMxyz,~,~]=Perform_Bspline_Transform(IMxyz,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
    end
    IM =  zeros(size(IM),'uint8');
    IMind=sub2ind(size(IM),round(IMxyz(1,:)),round(IMxyz(2,:)),round(IMxyz(3,:)));
    IM(IMind)=Int;
    IM = imtranslate(IM,[0, TR]);
    IM=max(IM,[],3);
    IMAll = max(IMAll,IM);
    figure(N);imshow(IMAll,[0 20])
    drawnow
    TR = TR - AddSpace;
%     
%     for j=ID-1:-1:1
%         StackPosition_prime=[1;1;1];
%         [IM,StackPosition_prime,~]=Perform_Bspline_Transform(IM,StackPosition_prime,L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
%     end
%     
    
%     IM_NR_max=max(IM,[],3);
%     IM_NR_max=imtranslate(IM_NR_max,StackPosition_prime(1:2));
%     IM_NR_max=IM_NR_max(1:1024,1:1024);
%     
% %     MIN=min([1;1],StackPosition_prime(1:2));
% %     MAX=max([1024;1024],size(IM_target_max)'+StackPosition_prime(1:2)-1);
% %     IMmax=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');
% %     IMmax(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_NR_max,1),...
% %     StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_NR_max,2))=IM_NR_max;
%     IM_NR_max = imtranslate(IM_NR_max,[0, TR]);
%     %             IMmax = imresize(IMmax,[1024,1024]);
%     IMAll = max(IMAll,IM_NR_max);
%     figure(N);imshow(IMAll,[0 20])
%     drawnow
%     TR = TR - AddSpace;
end
figure(N);imshow(IMAll,[0 20])
save('C:\Users\Seyed\Documents\Meetings\Research\SPIE\Results\AllFigs_trace2.mat','IMAll');








for ID = 1: size(StackList,1)
    sourceID = ID;
    targetID = ID + 1;
    
    Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';
    [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
    
    fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
    fname={fname.name}';
    Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
    [AM_source,r_source,R]=swc2AM(Path);
    [AM_source,r_source,~] = AdjustPPM(AM_source,r_source,R,ppm);
    
    
    if ID < size(StackList,1)
        
        fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
        fname={fname.name}';
        Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
        [AM_target,r_target,R]=swc2AM(Path);
        [AM_target,r_target,~] = AdjustPPM(AM_target,r_target,R,ppm);
        
        B1 = Boutons{sourceID,1};
        SourcePoints = B1.r1;
        TargetPoints = B1.r2;
        
        if ID > 1
            %         sourcePoints_NR_temp = SourcePoints';
            %         r_source_NR_temp = r_source';
            %         for j=size(b,2):-1:1
            %             [sourcePoints_NR_temp,~]=Perform_Bspline_Transform(sourcePoints_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
            %             [r_source_NR_temp,~]=Perform_Bspline_Transform(r_source_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
            %         end
            %         SourcePoints = sourcePoints_NR_temp';
            %         r_source = r_source_NR_temp';
            r_source=[];SourcePoints=[];AM_source=[];
            r_source = r_targetAll{sourceID-1};
            SourcePoints = TargetPointsAll{sourceID-1};
            AM_source = AM_targetAll{sourceID-1};
        end
        
        targetPoints_NR_temp = TargetPoints';
        r_target_NR_temp = r_target';
        for j=size(b,2):-1:1
            [targetPoints_NR_temp,~]=Perform_Bspline_Transform(targetPoints_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
            [r_target_NR_temp,~]=Perform_Bspline_Transform(r_target_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
        end
        TargetPoints = targetPoints_NR_temp';
        r_target = r_target_NR_temp';
        r_targetAll{sourceID} = r_target;
        AM_targetAll{sourceID}=AM_target;
        TargetPointsAll{sourceID} = TargetPoints;
        
        [k,f] = dsearchn(r_source,SourcePoints);
        
        s = 1;
        SourcePoints_onTrace=[];
        for j=1:size(k)
            if f(j)< 3
                SourcePoints_onTrace(s,:) = SourcePoints(j,:);
                s = s +1;
            end
        end
        figure(N); hold on;
        %     if ID==1
        %     figure(N); hold on; PlotAM(AM_source,r_source,'c')
        %     end
        
        plot3(SourcePoints_onTrace(:,2),SourcePoints_onTrace(:,1)+TR,ones(1,length(SourcePoints_onTrace(:,2))),'or')
        
        TR = TR - AddSpace;
        [k,f] = dsearchn(r_target,TargetPoints);
        TargetPoints(:,1) = TargetPoints(:,1)+TR;
        r_target(:,1)=r_target(:,1)+TR;
        %     PlotAM(AM_target,r_target,'c')
        s = 1;
        TargetPoints_onTrace=[];
        for j=1:size(k)
            if f(j)< 3
                TargetPoints_onTrace(s,:) = TargetPoints(j,:);
                s = s +1;
            end
        end
        plot3(TargetPoints_onTrace(:,2),TargetPoints_onTrace(:,1),ones(1,length(TargetPoints_onTrace(:,2))),'or')
        
        %         if ID==1
        %             line([SourcePoints_onTrace(:,2) TargetPoints_onTrace(:,2)]', ...
        %                 [SourcePoints_onTrace(:,1) TargetPoints_onTrace(:,1)]',[1;1]*ones(1,length(TargetPoints_onTrace(:,1))), 'Color', 'c');
        %         else
        %             line([SourcePoints_onTrace(:,2) TargetPoints_onTrace(:,2)]', ...
        %                 [SourcePoints_onTrace(:,1)+TR+AddSpace TargetPoints_onTrace(:,1)]',[1;1]*ones(1,length(TargetPoints_onTrace(:,1))), 'Color', 'c');
        %         end
        SourcePoints_onTrace;
        TargetPoints_onTrace;
        HC = [];
        for i=1:size(SourcePoints_onTrace,1)
            for j=1:size(TargetPoints_onTrace,1)
                HC(i,j) = mean(sum((SourcePoints_onTrace(i,:)-TargetPoints_onTrace(j,:)).^2,2).^0.5);
            end
        end
        
        Am = Hungarian_fast(HC);
        [idx1,idx2]=find(Am);
        x_Source = SourcePoints_onTrace(idx1,1);
        x_Target = TargetPoints_onTrace(idx2,1);
        y_Source = SourcePoints_onTrace(idx1,2);
        y_Target = TargetPoints_onTrace(idx2,2);
        z_Source = SourcePoints_onTrace(idx1,3);
        z_Target = TargetPoints_onTrace(idx2,3);
        matchLoc_Source = [x_Target,y_Target,z_Target];
        matchLoc_Target = [x_Source,y_Source,z_Source];
        
        if ID==1
        line([matchLoc_Source(:,2) matchLoc_Target(:,2)]', ...
                [matchLoc_Source(:,1) matchLoc_Target(:,1)]',...
                [1;1]*ones(1,length(matchLoc_Source(:,1))),'Color', rand(1,3));
        else
         line([matchLoc_Source(:,2) matchLoc_Target(:,2)]', ...
                [matchLoc_Source(:,1) matchLoc_Target(:,1)+TR+AddSpace]',...
                [1;1]*ones(1,length(matchLoc_Source(:,1))),'Color', rand(1,3));  
        end
        axis equal
        
        
        
        
        
        
        
        drawnow
    end
    
end










%%
% openfig('E:\Fig4B.fig');
L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
TR = -600;
for sourceID = 1: size(StackList,1)-1

    targetID = sourceID + 1;
    Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';
    [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
    [targetPoints_NR_temp,~]=Perform_Bspline_Transform(targetPoints',[],L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID},affine);
    
    
    
    B1 = Boutons{sourceID,1};
    SourcePoints = B1.r1;
    targetPoints = B1.r2;
    
%     for j=size(b,2):-1:1
%     [targetPoints_NR_temp,~]=Perform_Bspline_Transform(targetPoints',[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
%     end
    targetPoints_NR_temp = targetPoints_NR_temp';
    
    fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
    fname={fname.name}';
    Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
    [AM,r,R]=swc2AM(Path);
    [AM,r,~] = AdjustPPM(AM,r,R,ppm);
    [k,f] = dsearchn(r,SourcePoints);
    s = 1;
    for j=1:size(k)
        if f(j)< 3
%             figure(N);hold on; plot(SourcePoints(j,2)+200,SourcePoints(j,1),'or')
            SourcePoints_onTrace(s,:) = SourcePoints(j,:);
            s = s +1;
        end
    end
    
    fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
    fname={fname.name}';
    Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
    [AM,r,R]=swc2AM(Path);
    [AM,r,~] = AdjustPPM(AM,r,R,ppm);
    [k_target,f_target] = dsearchn(r,targetPoints_NR_temp);
    s = 1;
    for j=1:size(k)
        if f(j)< 3
%             figure(N);hold on; plot(SourcePoints(j,2)+200,SourcePoints(j,1),'or')
            TargetPoints_NR_onTrace(s,:) = targetPoints_NR_temp(j,:);
            s = s +1;
        end
    end

    
%     hold on; PlotAM(AM,Points_NR_temp',color(TraceNum+sourceID,:))
    
    for i=1:size(SourcePoints_onTrace,1)
        for j=1:size(TargetPoints_NR_onTrace,1)
            HC(i,j) = mean(sum((SourcePoints_onTrace(i,:)-TargetPoints_NR_onTrace(j,:)).^2,2).^0.5);
        end
    end
    TR = TR + 50;
    
    Am = Hungarian_fast(HC);
    [idx1,idx2]=find(Am);
    x_Source = SourcePoints_onTrace(idx1,1)+TR;
    x_Target = TargetPoints_NR_onTrace(idx2,1);
    y_Source = SourcePoints_onTrace(idx1,2)+200;
    y_Target = TargetPoints_NR_onTrace(idx2,2);
    z_Source = SourcePoints_onTrace(idx1,3);
    z_Target = TargetPoints_NR_onTrace(idx2,3);
    matchLoc_Target = [x_Target,y_Target,z_Target];
    matchLoc_Source = [x_Source,y_Source,z_Source];
    for i=1:size(matchLoc_Source,1)
        figure(N);hold on;line([matchLoc_Source(i,2) matchLoc_Target(i,2)], ...
            [matchLoc_Source(i,1)+200 matchLoc_Target(i,1)+200], 'Color', rand(1,3));
        
    end
    axis equal


end









for ID = 1: size(StackList,1)
    sourceID = ID - 1;
    targetID = ID 
    
    
    
    fname = dir([GTpath,'Matches\Traces\DL083',T_Names{ID},'001-A0*']);
    fname={fname.name}';
    Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
    [AM,r,R]=swc2AM(Path);
    [AM,r,~] = AdjustPPM(AM,r,R,ppm);
    Stack_File = char(StackList(ID,1));
    IM=ImportStack(char(Stack_File));
    IM = uint8(double(IM)./double(max(IM(:))).*255);
    IM_max=max(IM,[],3);
    
    
    for i=1:size(SourcePoints_onTrace,1)
        for j=1:size(TargetPoints_NR_onTrace,1)
            HC(i,j) = mean(sum((SourcePoints_onTrace(i,:)-TargetPoints_NR_onTrace(j,:)).^2,2).^0.5);
        end
    end
    Am = Hungarian_fast(HC);
    [idx1,idx2]=find(Am);
    x_Source = SourcePoints_onTrace(idx1,1);
    x_Target = TargetPoints_NR_onTrace(idx2,1);
    y_Source = SourcePoints_onTrace(idx1,2);
    y_Target = TargetPoints_NR_onTrace(idx2,2);
    z_Source = SourcePoints_onTrace(idx1,3);
    z_Target = TargetPoints_NR_onTrace(idx2,3);
    matchLoc_Target = [x_Target,y_Target,z_Target];
    matchLoc_Source = [x_Source,y_Source,z_Source];
    for i=1:size(matchLoc_Source,1)
        figure(N);hold on;line([matchLoc_Source(i,2)+200 matchLoc_Target(i,2)+200], ...
            [matchLoc_Source(i,1) matchLoc_Target(i,1)], 'Color', rand(1,3));
        
    end
    axis equal
    
end












%                     ----------------- Fig 4.B
[~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);

fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
fname={fname.name}';
Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
[AM,r,R]=swc2AM(Path);
[AM,r,~] = AdjustPPM(AM,r,R,ppm);
figure(N); hold on; PlotAM(AM,r,color(TraceNum+sourceID-1,:))
[k,f] = dsearchn(r,SourcePoints);
s = 1;
for j=1:size(k)
    if f(j)< 3
        figure(N);hold on; plot(SourcePoints(j,2),SourcePoints(j,1),'ob')
        SourcePoints_onTrace(s,:) = SourcePoints(j,:);
        s = s +1;
    end
end

fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
fname={fname.name}';
Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
[AM,r,R]=swc2AM(Path);
[AM,r,~] = AdjustPPM(AM,r,R,ppm);
[r_NR,~]=Perform_Bspline_Transform(r',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
r_NR = r_NR';
figure(N); hold on; PlotAM(AM,r_NR,color(TraceNum+sourceID-1,:))
[TargetPoints_NR,~]=Perform_Bspline_Transform(TargetPoints',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
TargetPoints_NR = TargetPoints_NR';
%     mean(mean((BTargetPoints-r_NR).^2,2).^0.5)
[k,f] = dsearchn(r_NR,TargetPoints_NR);
s = 1;
for j=1:size(k)
    if f(j)< 3
        figure(N);hold on; plot(TargetPoints_NR(j,2),TargetPoints_NR(j,1),'+g')
        TargetPoints_NR_onTrace(s,:) = TargetPoints_NR(j,:);
        s = s +1;
    end
end
% - before Registration
figure(N); hold on; PlotAM(AM,r,color(TraceNum+targetID,:))
for j=1:size(k)
    if f(j)< 3
        figure(N);hold on; plot(TargetPoints(j,2),TargetPoints(j,1),'+r')
    end
end

for i=1:size(SourcePoints_onTrace,1)
    for j=1:size(TargetPoints_NR_onTrace,1)
    HC(i,j) = mean(sum((SourcePoints_onTrace(i,:)-TargetPoints_NR_onTrace(j,:)).^2,2).^0.5);
    end
end
Am = Hungarian_fast(HC);
[idx1,idx2]=find(Am);
x_Source = SourcePoints_onTrace(idx1,1);
x_Target = TargetPoints_NR_onTrace(idx2,1);
y_Source = SourcePoints_onTrace(idx1,2);
y_Target = TargetPoints_NR_onTrace(idx2,2);
z_Source = SourcePoints_onTrace(idx1,3);
z_Target = TargetPoints_NR_onTrace(idx2,3);
matchLoc_Target = [x_Target,y_Target,z_Target];
matchLoc_Source = [x_Source,y_Source,z_Source];
for i=1:size(matchLoc_Source,1)
figure(N);hold on;line([matchLoc_Source(i,2) matchLoc_Target(i,2)], ...
                    [matchLoc_Source(i,1) matchLoc_Target(i,1)], 'Color', rand(1,3));    

end
axis equal

%% --------------------------------------------------- Hungarian Tracking
fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
fname={fname.name}';
Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
[AM,r,R]=swc2AM(Path);
[AM,r,~] = AdjustPPM(AM,r,R,ppm);
% figure(N); hold on; PlotAM(AM,r,color(TraceNum+sourceID-1,:))
[k,f] = dsearchn(r,SourcePoints);
for j=1:size(k)
%     if f(j)< 3
        figure(N);hold on; plot(SourcePoints(j,2),SourcePoints(j,1),'or')
%     end
end

fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
fname={fname.name}';
Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
[AM,r,R]=swc2AM(Path);
[AM,r,~] = AdjustPPM(AM,r,R,ppm);
[r_NR,~]=Perform_Bspline_Transform(r',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
r_NR = r_NR';
% figure(N); hold on; PlotAM(AM,r_NR,color(TraceNum+sourceID-1,:))
[TargetPoints_NR,~]=Perform_Bspline_Transform(TargetPoints',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
TargetPoints_NR = TargetPoints_NR';
%     mean(mean((BTargetPoints-r_NR).^2,2).^0.5)
    
    [k,f] = dsearchn(r_NR,TargetPoints_NR);
    for j=1:size(k)
%         if f(j)< 3
            figure(N);hold on; plot(TargetPoints_NR(j,2),TargetPoints_NR(j,1),'+g')
%         end
    end
for i=1:size(SourcePoints,1)
    for j=1:size(TargetPoints_NR,1)
    HC(i,j) = mean(sum((SourcePoints(i,:)-TargetPoints_NR(j,:)).^2,2).^0.5);
    end
end
Am = Hungarian_fast(HC);
[idx1,idx2]=find(Am);
x_Source = SourcePoints(idx1,1);
x_Target = TargetPoints_NR(idx2,1);
y_Source = SourcePoints(idx1,2);
y_Target = TargetPoints_NR(idx2,2);
z_Source = SourcePoints(idx1,3);
z_Target = TargetPoints_NR(idx2,3);
matchLoc_Target = [x_Target,y_Target,z_Target];
matchLoc_Source = [x_Source,y_Source,z_Source];
for i=1:size(matchLoc_Source,1)
figure(N);hold on;line([matchLoc_Source(i,2) matchLoc_Target(i,2)], ...
                    [matchLoc_Source(i,1) matchLoc_Target(i,1)], 'Color', rand(1,3));
result(i) = mean(sum((matchLoc_Source(i,:)-matchLoc_Target(i,:)).^2,2).^0.5);      
if result(i)<1
   x(i,1:2) = [result(i),1];
else
    x(i,1:2) = [result(i),0];
end
end
% TP =  size(find(result<=3),2);
% FP =  size(find(result>3),2);
roc(x)





































%%                     ----------------- Fig 4.A OLD
% IMAll = uint8(zeros(1024,1024));
% TR = -600;
% for i = 1: size(StackList,1)
%     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{i},'001-A0*']);
%     fname={fname.name}';
%     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%     [AM,r,R]=swc2AM(Path);
%     [AM,r,~] = AdjustPPM(AM,r,R,ppm);
%     Stack_File = char(StackList(i,1));
%     IM=ImportStack(char(Stack_File));
%     IM = uint8(double(IM)./double(max(IM(:))).*255);
%     IM_max=max(IM,[],3);
%     
%     [KT]=FastMarchingTube(size(IM),r,3,[1,1,1]);
%     IMmax =  max(uint8(KT).*IM,[],3);
%     IMmax = imtranslate(IMmax,[200, 0]);
%     r(:,2) = r(:,2) + 200;
%     
%     TR = TR + 50;
%     IMmax = imtranslate(IMmax,[0, TR]);
%     r(:,1) = r(:,1) + TR;
%     
%     IMAll = max(IMAll,IMmax);
%     figure(N);imshow(IMAll,[0 20])
%     
% %     figure(N); hold on; PlotAM(AM,r,color(TraceNum+i-1,:))
% %    
% %     if i == size(StackList,1)
% %         B1 = Boutons{i-1,1};
% %         SourcePoints = B1.r2;
% %         
% %     else
% %         B1 = Boutons{i,1};
% %         SourcePoints = B1.r1;
% %         
% %     end
% %     
% %     SourcePoints(:,2)+TR;
% %     
% %     [k,f] = dsearchn(r,SourcePoints);
% %     for j=1:size(k)
% %         if f(j)< 3
% %             figure(N);hold on; plot(SourcePoints(j,2),SourcePoints(j,1)+200,'or')
% %         end
% %     end
% end




