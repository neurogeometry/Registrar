clear all
clc
close all
% StackList_csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\slices2.csv';
 StackList_csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\slicesHoltmatMess.csv';
%StackList_csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\slicesHoltmat.csv';
addpath('C:\Users\Seyed\Documents\DatasetTests\registrar\registrar\Functions');
StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));
[PathStr,FolderName]=fileparts(StackList_csv_pth);
% DataFolder=[PathStr,'/Results-',FolderName];
DataFolder='data\Fig4';
Stack = [];
StackRegistered = [];
T_Translation = load([DataFolder,'\T_Translation.mat']);
T_Rigid = load([DataFolder,'\T_Rigid.mat']);
T_Affine = load([DataFolder,'\T_Affine.mat']);
% %                        Create Mess
% options.overwrite = 1;
% for i = 1:size(StackList,1)
%     IM = imread(char(StackList(i,1)));
%     [PathStr,Filename,ext]=fileparts(char(StackList(i,1)));
%     if mod(i,2) == 0
%         T = [2 -2];
%     else
%         T = [-2 2];
%     end
%         
%     IM = imtranslate(IM,T);
%     StackMess(:,:,i) = IM; 
% 
%     saveastiff(IM, ['E:\SliceRegistration\DL083B001G_Mess\',Filename,'.tif'],options);
% end
% figure,imshow(max(StackMess,[],3),[0 100]);

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
% 
% % Rigid
% for i = 1:size(StackList,1)
%     IM = imread(char(StackList(i,1)));
% %     Stack(:,:,i) = IM; 
%     [IMs{i},StackPosition_registered(i,:)]=Perform_Linear_Transform(IM,[1,1,1],T_Rigid.T.L(:,3*(i-1)+1:3*(i-1)+3),T_Rigid.T.b(:,i));
%     Sizes(i,:)=[size(IMs{i}),1];
% end
% StackPosition_registered=round(StackPosition_registered);
% Start=min(StackPosition_registered);
% End=max(StackPosition_registered+Sizes-1);
% Shift=1-Start;
% StackRegistered_Rigid=zeros(End-Shift,class(IM));
% for i = 1:size(StackList,1)
%     StackRegistered_Rigid(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
%         StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
% end
% 
% 
% Affine
for i = 1:size(StackList,1)
    IM = imread(char(StackList(i,1)));
     Stack(:,:,i) = IM; 
    [IMs{i},StackPosition_registered(i,:)]=Perform_Linear_Transform(IM,[1,1,1],T_Affine.T.L(:,3*(i-1)+1:3*(i-1)+3),T_Affine.T.b(:,i));
    Sizes(i,:)=[size(IMs{i}),1];
end
StackPosition_registered=round(StackPosition_registered);
Start=min(StackPosition_registered);
End=max(StackPosition_registered+Sizes-1);
Shift=1-Start;
StackRegistered_Affine=zeros(End-Shift,class(IM));
for i = 1:size(StackList,1)
    StackRegistered_Affine(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
        StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
end

correl_Affine = [];
for i = 2:size(StackList,1)
    % Affine
    cor3 = corrcoef(double(StackRegistered_Affine(:,:,i)),double(StackRegistered_Affine(:,:,i-1)));
    correl_Affine(i-1) = cor3(1,2);
    
end
% figure,imshow(max(StackRegistered_Affine,[],3),[0 100]);title('Affine');
k = 0;
mult = 400;
for j = 90:size(StackRegistered_Affine,3)-1
    
    k = k + 1;
    figure(k),imshowpair(uint16(Stack(:,:,j+1)).*mult,uint16(Stack(:,:,j)).*mult,'Scaling','independent')
    k = k + 1;
    figure(k),imshowpair(StackRegistered_Affine(:,:,j).*mult,StackRegistered_Affine(:,:,j+1).*mult,'Scaling','independent')

end
% Affine
FeaturePositions_NR = load([DataFolder,'\MatchedPoints_Affine.mat']);%MatchedPoints_Non-Rigid.mat
L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
mu = 20;
affine = 0;
nxyz = [1024;1024;1];
Sizes = [];
for i = 1:size(StackList,1)
    
    
    IM = imread(char(StackList(i,1)));
%     Stack(:,:,i) = IM; 
%     [IMs{i},StackPosition_registered(i,:)]=Perform_Linear_Transform(IM,[1,1,1],T_Affine.T.L(:,3*(i-1)+1:3*(i-1)+3),T_Affine.T.b(:,i));
%     Sizes(i,:)=[size(IMs{i}),1];
    
    
    
    sourceID = i
    targetID = i + 1;

    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
    
    IMs{1} = imread(char(StackList(sourceID,1)));
    Sizes(:,1)=[size(IMs{1}),1];
    StackPosition_registered(:,1) = [1,1,1];
    IM_target = imread(char(StackList(targetID,1)));
    
    [IMs{2},StackPosition_registered(:,2),~]=Perform_Bspline_Transform(IM_target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    Sizes(:,2)=[size(IMs{2}),1];
    
    
    
    
    
    
    
    
end




StackPosition_registered=round(StackPosition_registered);
Start=min(StackPosition_registered);
End=max(StackPosition_registered+Sizes-1);
Shift=1-Start;
StackRegistered_Affine=zeros(End-Shift,class(IM));
for i = 1:size(StackList,1)
    StackRegistered_Affine(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
        StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
end

correl_Affine = [];
for i = 2:size(StackList,1)
    % Affine
    cor3 = corrcoef(double(StackRegistered_Affine(:,:,i)),double(StackRegistered_Affine(:,:,i-1)));
    correl_Affine(i-1) = cor3(1,2);
    
end



% Non-Rigid Pairwise works
FeaturePositions_NR = load([DataFolder,'\MatchedPoints_Affine.mat']);%MatchedPoints_Non-Rigid.mat
L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
mu = 20;
affine = 0;
nxyz = [1024;1024;1];
numIMGs = size(StackList,1)-1;
for i = 1:numIMGs
    sourceID = i
    targetID = i + 1;

    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
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
    StackRegistered_NR=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');

    for j = 1:2
    StackRegistered_NR(StackPosition_registered(j,1)-MIN(1)+1:StackPosition_registered(j,1)-MIN(1)+Sizes(1,j),...
                    StackPosition_registered(j,2)-MIN(2)+1:StackPosition_registered(j,2)-MIN(2)+Sizes(2,j),j)=IMs{j};
    end
%      figure,imshow(max(StackRegistered_NR,[],3),[0 max(StackRegistered_NR(:))]);
    cor3 = corrcoef(double(StackRegistered_NR(:,:,1)),double(StackRegistered_NR(:,:,2)));
    correl_NR(i) = cor3(1,2)
    
    
%     Start=min(StackPosition_registered);
%     End=max(StackPosition_registered+Sizes-1);
%     Shift=1-Start;
%     StackRegistered_NR=zeros(End-Shift,class(IMs{1}));
%     for j = 1:2
%     StackRegistered_NR(StackPosition_registered(j,1)+Shift(1):StackPosition_registered(j,1)+Shift(1)+Sizes(j,1)-1,...
%         StackPosition_registered(j,2)+Shift(2):StackPosition_registered(j,2)+Shift(2)+Sizes(j,2)-1,j)=IMs{j};
%     end
%     cor3 = corrcoef(double(StackRegistered_NR(:,:,1)),double(StackRegistered_NR(:,:,2)));
%     correl_NR1(i) = cor3(1,2)
    
   
end
    
% load('C:\Users\Seyed\Documents\Meetings\Research\result_withinStack\correl_NR.mat');
mean(correl_NR)
correl_NonRigid = correl_NR;




correl_before=[];
correl_Translation=[];
correl_Rigid=[];
correl_Affine =[];

for i = 2:size(StackList,1)
    cor3 = corrcoef(double(Stack(:,:,i)),double(Stack(:,:,i-1)));
    correl_before(i-1) = cor3(1,2);
    % Translation
    cor3 = corrcoef(double(StackRegistered_Translation(:,:,i)),double(StackRegistered_Translation(:,:,i-1)));
    correl_Translation(i-1) = cor3(1,2);
    % Rigid
    cor3 = corrcoef(double(StackRegistered_Rigid(:,:,i)),double(StackRegistered_Rigid(:,:,i-1)));
    correl_Rigid(i-1) = cor3(1,2);
    % Affine
    cor3 = corrcoef(double(StackRegistered_Affine(:,:,i)),double(StackRegistered_Affine(:,:,i-1)));
    correl_Affine(i-1) = cor3(1,2);
%     if i < 100
%     cor3 = corrcoef(double(StackRegistered_NR(:,:,i)),double(StackRegistered_NR(:,:,i-1)));
%     correl_NR(i-1) = cor3(1,2);
%     end
    
end

mean(correl_before)
std(correl_before)
mean(correl_Translation)
std(correl_Translation)
mean(correl_Rigid)
std(correl_Rigid)
mean(correl_Affine)
std(correl_Affine)
mean(correl_NonRigid)
std(correl_NonRigid)

figure,
boxplot([correl_before(:),correl_Translation(:),correl_Rigid(:),...
        correl_Affine(:),correl_NonRigid(:)],'Whisker',inf)
    axis square, box on

figure,
plot(correl_before)
hold on
plot(correl_Translation)
plot(correl_Rigid)
plot(correl_Affine)
plot(correl_NonRigid)

figure,
plot(sort(correl_before))
hold on
plot(sort(correl_Translation))
plot(sort(correl_Rigid))
plot(sort(correl_Affine))
plot(sort(correl_NonRigid))

% for i = 1:size(StackList,1)
%     IM = imread(char(StackList(i,1)));
%     Stack(:,:,i) = IM; 
%     IM = imtranslate(IM,T_Translation.T.b([2,1],i)');
%     StackRegistered(:,:,i) = IM; 
% end
% figure,imshow(max(Stack,[],3),[0 max(Stack(:))]);
% figure,imshow(max(StackRegistered,[],3),[0 max(StackRegistered(:))]);

figure,imshow(max(Stack,[],3),[0 100]);title('Original');
figure,imshow(max(StackRegistered_Translation,[],3),[0 100]);title('Translation');
figure,imshow(max(StackRegistered_Rigid,[],3),[0 100]);title('Rigid');
figure,imshow(max(StackRegistered_Affine,[],3),[0 100]);title('Affine');
% figure,imshow(max(StackRegistered_NR,[],3),[0 100]);title('NonRigid');













% -------------------------- Non-Rigid Works
% FeaturePositions_NR = load([DataFolder,'\MatchedPoints_Non-Rigid.mat']);
% L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
% IMs = [];
% TR = 0;
% AddSpace = 40;
% mu = 1;
% affine = 0;
% nxyz = [256;256;1];
% for ID = 1: size(StackList,1)-2
%     sourceID = ID;
%     targetID = ID + 1;
% 
%     Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
%     Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
%     [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
% 
%             Stack_File = char(StackList(sourceID,1));
%             IM=ImportStack(char(Stack_File));
%             IM = uint8(double(IM)./double(max(IM(:))).*255);
%             
%             if ID > 1
%                 for j=ID-1:-1:1
%                 [IM,StackPosition_prime,~]=Perform_Bspline_Transform(IM,[1;1;1],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
%                 end
%                 IMmax=IM;%max(IM,[],3);
% %                 MIN=min([1;1],StackPosition_prime(1:2));
% %                 MAX=max(size(IM1)',size(IM)'-StackPosition_prime(1:2)+1);
% %                 temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');
% %                 IMmax=temp;
% %                 IMmax(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_NR_max,1),...
% %                     StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_NR_max,2))=IM_NR_max;
% %                 IMmax = IM_NR_max;
%                 StackPosition_registered(ID,:) = StackPosition_prime;
%             else
%                 IMmax =  IM;%max(IM,[],3);
%                 StackPosition_registered(ID,:) = [1,1,1];
%             end
% %             IMmax = imtranslate(IMmax,[0, TR]);
%  %           IMAll(:,:,ID) = IMmax;%max(IMAll,IMmax);
%              IMAll{ID} = IMmax;%max(IMAll,IMmax);
%              Sizes(ID,:)=[size(IMAll{ID}),1];
% %             TR = TR - AddSpace;
% end

% load('C:\Users\Seyed\Documents\Meetings\Research\result_withinStack\stackRegisteredNR.mat','IMAll','Sizes','StackPosition_registered','StackRegistered_NR');
% StackPosition_registered=round(StackPosition_registered);
% Start=min(StackPosition_registered);
% End=max(StackPosition_registered+Sizes-1);
% Shift=1-Start;
% StackRegistered_NR=zeros(End-Shift,class(IM));
% for i = 1:size(StackList,1)-2
%     StackRegistered_NR(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
%         StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMAll{i};
% end
% correl_NonRigid = [];
% for i = 2:size(StackList,1)-2
%     % Non-Rigid
%     cor3 = corrcoef(double(StackRegistered_NR(:,:,i)),double(StackRegistered_NR(:,:,i-1)));
%     correl_NonRigid(i-1) = cor3(1,2);
%     
% end
% correl_NonRigid(99:100) = 1;
% mean(correl_NonRigid)




% % % ------------------------ Non-Rigid
% % 
% % FeaturePositions_NR = load([DataFolder,'\MatchedPoints_Affine.mat']);
% % L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
% % IMs = [];
% % TR = 0;
% % AddSpace = 40;
% % mu = 1040;
% % affine = 1;
% % nxyz = [256;256;1];
% % numIMGs = size(StackList,1)
% % for ID = 1: numIMGs
% %     sourceID = ID;
% %     targetID = ID + 1;
% % 
% %     Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
% %     Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
% %     [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
% %             Stack_File = char(StackList(sourceID,1));
% %             IM = imread(char(StackList(sourceID,1)));
% %             org{ID} = IM;
% %             if ID > 1
% %                 for j=size(b,2)-1:-1:1
% % %                  j = sourceID - 1;
% %                 [IM,StackPosition_prime,~]=Perform_Bspline_Transform(IM,[1;1;1],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
% % %                 j
% %                 end
% %                 IMs{ID} = IM;
% %                 StackPosition_registered(ID,:) = StackPosition_prime;
% %                 
% %            else
% %                 IMs{ID} = IM;
% %                 StackPosition_registered(ID,:) = [1,1,1];
% %             end
% %             Sizes(ID,:)=[size(IMs{ID}),1];
% %             ID
% % end
% % StackPosition_registered=round(StackPosition_registered);
% % Start=min(StackPosition_registered);
% % End=max(StackPosition_registered+Sizes-1);
% % Shift=1-Start;
% % StackRegistered_NR=zeros(End-Shift,class(IM));
% % for i = 1:numIMGs
% %     StackRegistered_NR(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
% %         StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
% % end
% % 
% % 
% % 
% % figure,imshow(max(StackRegistered_NR,[],3),[0 100])
% % correl_NR = [];
% % for i = 2:5
% %     % Non-Rigid
% %     cor3 = corrcoef(double(StackRegistered_NR(:,:,i)),double(StackRegistered_NR(:,:,i-1)));
% %     correl_NR(i-1) = cor3(1,2);
% %     
% % end
% % mean(correl_NR)
% % i
% 
% 
% 
% 
% % L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
% % mu = 1040;
% % affine = 1;
% % nxyz = [256;256;156];
% % % AffineNonRigid
% % StackPosition_registered = [];
% % StackRegistered_Non_rigid=[];
% % for i = 1:size(StackList,1)-1
% %     sourceID = i;
% %     targetID = i + 1;
% %     
% %     
% %     FeaturePositions_NR = load([DataFolder,'\MatchedPoints_Non-Rigid.mat']);
% %     Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
% %     Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
% %     [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
% %     IM_Source = imread(char(StackList(sourceID,1)));
% %     if i > 1
% %         [IMs{sourceID},StackPosition_registered(:,sourceID),~]=Perform_Bspline_Transform(IM_Source,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
% %     else
% %         IMs{sourceID} = IM_Source;
% %     end
% %     Sizes(sourceID,:)=[size(IMs{sourceID}),1];
% % end
% % StackPosition_registered(1,1,3) =1;
% % StackPosition_registered=round(StackPosition_registered');
% % Start=min(StackPosition_registered);
% % End=max(StackPosition_registered+Sizes);
% % Shift=1-Start;
% % StackRegistered_Non_rigid=zeros(End-Shift,class(IM_Source));
% % for i = 1:size(StackList,1)-1
% % %     StackRegistered_Non_rigid(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
% % %         StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
% %     StackRegistered_Non_rigid(:,:,i) = IMs{i};
% % end

