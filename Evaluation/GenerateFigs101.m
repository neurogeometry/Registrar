clear all;
close all;
clc;
N = 1;
sourceID = 1;
targetID = 2;
% for k = 1:17
%     sourceID = k;
%     targetID = k + 1;
TraceNum = 2;
mu = 1040;
affine = 1;
ppm = 1;
 %Image Numbers
addpath('NeuronTracerV20');
csvPath = '..\..\RegistrationEvaluation\TimeLapse_Holtmaat_101_StackList.csv';
% for matching, comment otherwise
 csvPath = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\TimeLapse_Holtmaat_101_StackList_forMatchingResults.csv';

StackList = table2cell(readtable(csvPath,'Delimiter',','));
GTpath = 'E:\Datasets\TimeLaps\Mouse101\original\Profile\';
Functionspath = '..\';

resultpath = '..\..\RegistrationEvaluation\Results-TimeLapse_Holtmaat_101_StackList\';
% for matching , comment otherwise
 resultpath = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-TimeLapse_Holtmaat_101_StackList_forMatchingResults\';

addpath([Functionspath,'Functions']);

% Laod all matched Boutons
% load([GTpath,'Matches\DL083-001-Matches.mat']);


T_Names = {'B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S'};
T_Days = 1:4:41;

% calculate Nonrigid Transform
%  FeaturePositions_NR = load([resultpath,'MatchedPoints_Non-Rigid_mu0.mat']);

%FeaturePositions_NR = load('MatchedPoints_Non-Rigid_mu1020.mat');
% Matched_hung = FeaturePositions_NR.Matched_hung;

% FeaturePositions_NR = load('..\..\RegistrationEvaluation\Results-TimeLapse_Holtmaat_101_StackList\MatchedPoints_Translation.mat');
FeaturePositions_NR = load([resultpath,'MatchedPoints_Non-Rigid.mat']);
Matched_hung = FeaturePositions_NR.Matched_hung;

Matched = FeaturePositions_NR.Matched;

% Matched_hung = FeaturePositions_NR.Matched_hung;
Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';

nxyz = [256;256;156];


%Use Boutons
% Boutons = struct2cell(Matches);
% B1 = Boutons{sourceID,1};
% SourcePoints = B1.r1;
% TargetPoints = B1.r2;
% BoutonsBefore = mean(sum((TargetPoints-SourcePoints).^2,2).^0.5)
% BoutonErrorBefore = mean(sum((TargetPoints-SourcePoints).^2,2).^0.5)/sqrt(size(TargetPoints,1))


% k=0;
% for i = 1:12
% B1 = Boutons{i,1};
% k = k + size(B1.r1,1);
% end

% Get all Traces
fname_First = dir([GTpath,'DL101',T_Names{sourceID},'008\']);
fname_First = fname_First(3:end);
fname_First=cell2mat({fname_First.name}');

fname_Second = dir([GTpath,'DL101',T_Names{targetID},'008\']);
fname_Second = fname_Second(3:end);
fname_Second=cell2mat({fname_Second.name}');
color = rand(size(fname_First,1)+50,3);

% Laod images
Source_Stack_File = char(StackList(sourceID,1));
Target_Stack_File = char(StackList(targetID,1));
IM_Source=ImportStack(char(Source_Stack_File));
IM_Source = uint8(double(IM_Source)./double(max(IM_Source(:))).*255);
IM_Target=ImportStack(char(Target_Stack_File));
IM_Target = uint8(double(IM_Target)./double(max(IM_Target(:))).*255);
IM_source_max=max(IM_Source,[],3);
IM_target_max=max(IM_Target,[],3);

% %% Show All Traces on Image  ----------------- Fig 2. A,B
% N = N + 1;
% figure(N);imshow(IM_source_max,[0 50]);
% N = N + 1;
% figure(N);imshow(IM_target_max,[0 50]);
% for i=1:size(fname_First,1) % All Traces
%     % Using Trace
% 
%     sourcePath = [GTpath,'DL101',T_Names{sourceID},'008\',fname_First(i,:)];
%     targetPath = [GTpath,'DL101',T_Names{targetID},'008\',fname_Second(i,:)];
% 
%     load(sourcePath)
%     AM_Source = AM.optim;
%     r_Source = r.optim;
%     R_Source = ones(size(r_Source,1),1);
%     clear 'AM'
%     clear 'r'
%     
%     load(targetPath)
%     AM_Target = AM.optim;
%     r_Target = r.optim;
%     R_Target = ones(size(r_Target,1),1);
% 
%     [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
%     [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
% 
%     figure(N-1);hold on; PlotAM(AM_Source,r_Source,color(i,:))
% %     text(r_Source(round(size(r_Source,1)/3),2),r_Source(round(size(r_Source,1)/3),1),['Axon #',num2str(i)],'Color',color(i,:))
%     figure(N);hold on; PlotAM(AM_Target,r_Target,color(i,:))
% %     text(r_Target(round(size(r_Target,1)/3),2),r_Target(round(size(r_Target,1)/3),1),['Axon #',num2str(i)],'Color',color(i,:))
% 
% end

%% Show One Trace  ----------------- Fig 2. A,B - Trace
% TraceNum = 6;
% % for TraceNum=1:22
% sourcePath = [GTpath,'DL101',T_Names{sourceID},'008\',fname_First(TraceNum,:)];
% targetPath = [GTpath,'DL101',T_Names{targetID},'008\',fname_Second(TraceNum,:)];
% load(sourcePath)
% AM_Source = AM.optim;
% r_Source = r.optim;
% R_Source = ones(size(r_Source,1),1);
% clear 'AM'
% clear 'r'
% load(targetPath)
% AM_Target = AM.optim;
% r_Target = r.optim;
% R_Target = ones(size(r_Target,1),1);
% 
% [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
% [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
% 
% [KT]=FastMarchingTube(size(IM_Source),r_Source,10,[1,1,1]);
% N = N+1;
% figure(N);imshow(max(uint8(KT).*IM_Source,[],3),[0 20])
% figure(N);hold on; PlotAM(AM_Source,r_Source,color(TraceNum,:))
% % figure(N);hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'or')
% % text(r_Source(round(size(r_Source,1)/3),2),r_Source(round(size(r_Source,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(0)],'Color',color(TraceNum,:))
% 
% [KT]=FastMarchingTube(size(IM_Target),r_Target,10,[1,1,1]);
% N = N+1;
% figure(N);imshow(max(uint8(KT).*IM_Target,[],3),[0 20])
% figure(N);hold on; PlotAM(AM_Target,r_Target,color(TraceNum,:))
% figure(N);hold on; plot(TargetPoints(:,2),TargetPoints(:,1),'or')
% text(r_Target(round(size(r_Target,1)/3),2),r_Target(round(size(r_Target,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(4)],'Color',color(TraceNum,:))
% end

%% Show Overlaped Images Fig. 2. C
% N = N+1;
% figure(N),imshowpair(IM_source_max,IM_target_max,'Scaling','independent')

% Show One Trace All sessions  ----------------- Fig 2. D
%
% for TraceNum = 12:12
%      TraceNum 
% for i = 1: size(StackList,1)
%     
%     
%     sourcePath = [GTpath,'DL101',T_Names{i},'008\',fname_First(TraceNum,:)];
%     load(sourcePath)
%     AM = AM.optim;
%     r = r.optim;
%     R = ones(size(r,1),1);
% 
% %     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{i},'001-A0*']);
% %     fname={fname.name}';
% %     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
% %     [AM,r,R]=swc2AM(Path);
%     
%     
%     [AM,r,~] = AdjustPPM(AM,r,R,ppm);
%     if i ==1
%     Stack_File = char(StackList(i,1));
%     IM=ImportStack(char(Stack_File));
%     IM = uint8(double(IM)./double(max(IM(:))).*255);
%     IM_max=max(IM,[],3);
% 
%     [KT]=FastMarchingTube(size(IM),r,3,[1,1,1]);
%      N = N+1;
%     figure(N);imshow(max(uint8(KT).*IM,[],3),[0 20])
%     end
%     figure(N); hold on; PlotAM(AM,r,color(TraceNum+i-1,:))
% %     text(r(round(size(r,1)/3),2),r(round(size(r,1)/3),1),['Day ',num2str(T_Days(i))],'Color',color(TraceNum+i-1,:))
% %     text(r(round(size(r,1)/3),2),r(round(size(r,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(T_Days(i))],'Color',color(TraceNum,:))
% %     if i == size(StackList,1)
% %         B1 = Boutons{i-1,1};
% %         SourcePoints = B1.r2;
% % 
% %     else
% %         B1 = Boutons{i,1};
% %         SourcePoints = B1.r1;
% % 
% %     end
% %     figure(N);hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'or')
% end
% % saveas(figure(N),['C:\Users\Seyed\Documents\Meetings\Research\SPIE\Conference\Figs\Traces\',num2str(TraceNum),'_Before.fig']);
% end




% Matching Result  ----------------- Fig 2.
stitched = appendimages(IM_target_max,IM_source_max,'horizontal');
%                   ----------------- Fig 2. A
N = N + 1;
matchesStep = 6;
figure(N);imshow(stitched,[0 40])
figure(N);hold on;plot(Matched_hung{sourceID,targetID}(1:matchesStep:end,5)+1,Matched_hung{sourceID,targetID}(1:matchesStep:end,4)+1,'*r');
figure(N);hold on;plot(Matched_hung{sourceID,targetID}(1:matchesStep:end,2)+size(IM_target_max,1)+45+1,Matched_hung{sourceID,targetID}(1:matchesStep:end,1),'*r');
figure(N);hold on;
for i = 1:matchesStep:length(Matched_hung{sourceID,targetID})
line([Matched_hung{sourceID,targetID}(i,5)+1 Matched_hung{sourceID,targetID}(i,2)+size(IM_target_max,1)+45+1], ...
                    [Matched_hung{sourceID,targetID}(i,4)+1 Matched_hung{sourceID,targetID}(i,1)], 'Color', rand(1,3));
end

%                   ----------------- Fig 2. B
N = N + 1;
figure(N);imshow(stitched,[0 20])
figure(N);hold on;plot(Matched{sourceID,targetID}(:,5),Matched{sourceID,targetID}(:,4),'*r');
figure(N);hold on;plot(Matched{sourceID,targetID}(:,2)+size(IM_target_max,1)+45,Matched{sourceID,targetID}(:,1),'*r');
figure(N);hold on;
for i = 1: length(Matched{sourceID,targetID})
line([Matched{sourceID,targetID}(i,5) Matched{sourceID,targetID}(i,2)+size(IM_target_max,1)+45], ...
                    [Matched{sourceID,targetID}(i,4) Matched{sourceID,targetID}(i,1)], 'Color', rand(1,3));
end
i


%%                   ----------------- Fig 2. C Show Grid
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
[~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);

[IM_Target_NR,StackPosition_prime,~]=Perform_Bspline_Transform(IM_Target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
IM_Target_NR_max=max(IM_Target_NR,[],3);
MIN=min([1;1],StackPosition_prime(1:2));
MAX=max(size(IM_source_max)',size(IM_target_max)'+StackPosition_prime(1:2));
temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');

IM_Target_NR_max_P=temp;
IM_Target_NR_max_P(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_Target_NR_max,1),...
    StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_Target_NR_max,2))=IM_Target_NR_max;
IM_source_max_P=temp;
IM_source_max_P(2-MIN(1):1-MIN(1)+size(IM_source_max,1),...
    2-MIN(2):1-MIN(2)+size(IM_source_max,2))=IM_source_max;
N = N +1;
figure(N),imshowpair(IM_source_max_P,IM_Target_NR_max_P,'Scaling','independent')
% end
N
%%                    All Traces (From All Sessions) After Transform  ----------------- Fig 3.B
% 
% 
% for TraceNum = 4:20
% %      TraceNum 
%     L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
%     for sourceID = 1: size(StackList,1)-1
% 
% 
%         sourcePath = [GTpath,'DL101',T_Names{sourceID},'008\',fname_First(TraceNum,:)];
%         load(sourcePath)
%         AM = AM.optim;
%         r = r.optim;
%         R = ones(size(r,1),1);
% 
% 
%     %     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
%     %     fname={fname.name}';
%     %     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%     %     [AM,r,R]=swc2AM(Path);
%         [AM,r,~] = AdjustPPM(AM,r,R,ppm);
% 
%          if sourceID ==1
%             Stack_File = char(StackList(sourceID,1));
%         IM=ImportStack(char(Stack_File));
%         IM = uint8(double(IM)./double(max(IM(:))).*255);
%         IM_max=max(IM,[],3);
%         [KT]=FastMarchingTube(size(IM),r,3,[1,1,1]);
%         N = N+1;
%         figure(N);imshow(max(uint8(KT).*IM,[],3),[0 20])
%         figure(N);hold on; PlotAM(AM,r,color(TraceNum,:))
%     %     text(r(round(size(r,1)/3),2),r(round(size(r,1)/3),1),['Axon #',num2str(TraceNum),' Day ',num2str(T_Days(i))],'Color',color(TraceNum,:))
%     %     if i == size(StackList,1)
%     %         B1 = Boutons{i-1,1};
%     %         SourcePoints = B1.r2;
%     %     else
%     %         B1 = Boutons{i,1};
%     %         SourcePoints = B1.r1;
%     %     end
%     %     figure(N);hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'or')
%          end
% 
%         targetID = sourceID + 1;
%         Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
%         Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';
%         [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
% 
% 
% 
%         sourcePath = [GTpath,'DL101',T_Names{targetID},'008\',fname_First(TraceNum,:)];
%         load(sourcePath)
%         AM = AM.optim;
%         r = r.optim;
%         R = ones(size(r,1),1);
%     %     fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
%     %     fname={fname.name}';
%     %     Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
%     %     [AM,r,R]=swc2AM(Path);
% 
%         [AM,r,~] = AdjustPPM(AM,r,R,ppm);
%         Points_NR_temp = r';
%         for j=size(b,2):-1:1
%         [Points_NR_temp,~]=Perform_Bspline_Transform(Points_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
%         end
%         hold on; PlotAM(AM,Points_NR_temp',color(TraceNum+sourceID,:))
% 
%     %     [KT]=FastMarchingTube(size(IM),r,3,[1,1,1]);
%     %     Target_Stack_File = char(StackList(targetID,1));
%     %     IM_Target=ImportStack(char(Target_Stack_File));
%     %     IM_Target = uint8(double(IM_Target)./double(max(IM_Target(:))).*255);
%     %     IM_Target = uint8(KT).*IM_Target;
%     %     [IM_Target_NR,StackPosition_prime,~]=Perform_Bspline_Transform(IM_Target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
%     %     IM_Target_NR_max=max(IM_Target_NR,[],3);
%     %     MIN=min([1;1],StackPosition_prime(1:2));
%     %     MAX=max(size(IM_source_max)',size(IM_target_max)'+StackPosition_prime(1:2)-1);
%     %     temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');
%     %     IM_Target_NR_max_P=temp;
%     %     IM_Target_NR_max_P(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_Target_NR_max,1),...
%     %         StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_Target_NR_max,2))=IM_Target_NR_max;
%     %     IM_source_max_P=temp;
%     %     IM_source_max_P(2-MIN(1):1-MIN(1)+size(IM_source_max,1),...
%     %         2-MIN(2):1-MIN(2)+size(IM_source_max,2))=IM_source_max;
%     %     N = N +1;
%     %     figure(N),imshowpair(IM_Target_NR_max_P,IM_source_max_P,'Scaling','independent')
% 
%     end
% % end
% saveas(figure(N),['C:\Users\Seyed\Documents\Meetings\Research\SPIE\Conference\Figs\Traces\',num2str(TraceNum),'_After.fig']);
% end

%                      ----------------- Fig 3.C
%From Evaluate_TimeLapse.m file
%                      ----------------- Fig 3.D
%From Evaluate_TimeLapse.m file