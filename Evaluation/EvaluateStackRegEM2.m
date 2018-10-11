clear all;

StackList_csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\slice_DKB_7780_expl2_002-sub1.csv';
StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));
[PathStr,FolderName]=fileparts(StackList_csv_pth);
DataFolder=[PathStr,'/Results-',FolderName];

im1 = imread(char(StackList(1,1)));
im2 = imread(char(StackList(2,1)));
figure(1),imshowpair(im1,im2,'Scaling','independent')
cor3 = corrcoef(double(im1),double(im2));
correl_before(1) = cor3(1,2)


T_Translation = load([DataFolder,'\T_Translation.mat']);
T_Rigid = load([DataFolder,'\T_Rigid.mat']);
T_Affine = load([DataFolder,'\T_Affine.mat']);


% % Translation
% for i = 1:size(StackList,1)
%     IM = imread(char(StackList(i,1)));
%     Stack(:,:,i) = IM; 
%     [IMs{i},StackPosition_registered(:,i)]=Perform_Linear_Transform(IM,[1,1,1],[],T_Translation.T.b(:,i));
%     Sizes(i,:)=[size(IMs{i}),1];
% end
% StackPosition_registered=round(StackPosition_registered');
% Start=min(StackPosition_registered);
% End=max(StackPosition_registered+Sizes-1);
% Shift=1-Start;
% StackRegistered_Translation=zeros(End-Shift,class(IM));
% for i = 1:size(StackList,1)
%     StackRegistered_Translation(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
%         StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
% end
% figure(3),imshowpair(StackRegistered_Translation(:,:,1),StackRegistered_Translation(:,:,2),'Scaling','independent')
% figure(2),imshowpair(IMs{1},IMs{2},'Scaling','independent')

% Sizes = [];
% StackPosition_registered = [];
% % Rigid
% for i = 1:size(StackList,1)
%     IM = imread(char(StackList(i,1)));
% %     Stack(:,:,i) = IM; 
%     [IMs{i},StackPosition_registered(:,i)]=Perform_Linear_Transform(IM,[1,1,1],T_Rigid.T.L(:,3*(i-1)+1:3*(i-1)+3),T_Rigid.T.b(:,i));
%     Sizes(i,:)=[size(IMs{i}),1];
% end
% StackPosition_registered=round(StackPosition_registered');
% Start=min(StackPosition_registered);
% End=max(StackPosition_registered+Sizes-1);
% Shift=1-Start;
% StackRegistered_Rigid=zeros(End-Shift,class(IM));
% for i = 1:size(StackList,1)
%     StackRegistered_Rigid(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
%         StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
% end
% figure(3),imshowpair(StackRegistered_Rigid(:,:,1),StackRegistered_Rigid(:,:,2),'Scaling','independent')

% Sizes = [];
% StackPosition_registered = [];
% % Affine
% for i = 1:size(StackList,1)
%     if i == 1
%         j = 2;
%     else
%         j == 1;
%     end
%     IM = imread(char(StackList(j,1)));
% %     Stack(:,:,i) = IM; 
%     [IMs{i},StackPosition_registered(i,:)]=Perform_Linear_Transform(IM,[1,1,1],T_Affine.T.L(:,3*(i-1)+1:3*(i-1)+3),T_Affine.T.b(:,i));
%     Sizes(i,:)=[size(IMs{i}),1];
% end
% StackPosition_registered=round(StackPosition_registered);
% Start=min(StackPosition_registered);
% End=max(StackPosition_registered+Sizes-1);
% Shift=1-Start;
% StackRegistered_Affine=zeros(End-Shift,class(IM));
% for i = 1:size(StackList,1)
%     StackRegistered_Affine(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
%         StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
% end
% figure(4),imshowpair(StackRegistered_Affine(:,:,1),StackRegistered_Affine(:,:,2),'Scaling','independent')
% 
% cor3 = corrcoef(double(StackRegistered_Affine(:,:,1)),double(StackRegistered_Affine(:,:,2)));
% correl_Affine(1) = cor3(1,2)   

Sizes = [];
StackPosition_registered = [];
% Non-Rigid Pairwise works
FeaturePositions_NR = load([DataFolder,'\MatchedPoints_Affine.mat']);%MatchedPoints_Non-Rigid.mat
L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
mu = 800;
affine = 0;
nxyz = [64;64;1];
numIMGs = size(StackList,1)-1;
for i = 1:numIMGs
    sourceID = i
    targetID = i + 1;

    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
    
    IMs{1} = imread(char(StackList(sourceID,1)));
    Sizes(:,1)=[size(IMs{1}),1];
    StackPosition_registered(:,1) = [1,1,1];
    IM_target = imread(char(StackList(targetID,1)));
    
    [IMs{2},StackPosition_registered(:,2),~]=Perform_Bspline_Transform(IM_target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    Sizes(:,2)=[size(IMs{2}),1];
    
    StackPosition_registered=round(StackPosition_registered);

   
    StackRegistered_NR = [];
    MIN=min(min([1;1],StackPosition_registered(1:2,:)));
    MAX=max(max(size(IMs{1})',size(IMs{1})'+StackPosition_registered(1:2,:)-1));
    StackRegistered_NR=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint16');

    for j = 1:2
    StackRegistered_NR(StackPosition_registered(j,1)-MIN(1)+1:StackPosition_registered(j,1)-MIN(1)+Sizes(1,j),...
                    StackPosition_registered(j,2)-MIN(2)+1:StackPosition_registered(j,2)-MIN(2)+Sizes(2,j),j)=IMs{j};
    end
%      figure,imshow(max(StackRegistered_NR,[],3),[0 max(StackRegistered_NR(:))]);
    cor3 = corrcoef(double(StackRegistered_NR(:,:,1)),double(StackRegistered_NR(:,:,2)));
    correl_NR(i) = cor3(1,2)    
   
end

figure(5),imshowpair(StackRegistered_NR(:,:,1),StackRegistered_NR(:,:,2),'Scaling','independent')
% figure(6),imshow(max(StackRegistered_NR,[],3),[0 max(StackRegistered_NR(:))]);
i