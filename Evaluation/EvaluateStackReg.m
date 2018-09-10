clear all
clc
close all
% StackList_csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\slices2.csv';
 StackList_csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\slicesHoltmatMess.csv';
%StackList_csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\slicesHoltmat.csv';

StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));
[PathStr,FolderName]=fileparts(StackList_csv_pth);
DataFolder=[PathStr,'/Results-',FolderName];
Stack = [];
StackRegistered = [];
T_Translation = load([DataFolder,'\T_Translation.mat']);
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


% Translation
for i = 1:size(StackList,1)
    IM = imread(char(StackList(i,1)));
    Stack(:,:,i) = IM; 
    [IMs{i},StackPosition_registered(:,i)]=Perform_Linear_Transform(IM,[1,1,1],[],T_Translation.T.b(:,i));
    Sizes(i,:)=[size(IMs{i}),1];
end
StackPosition_registered=round(StackPosition_registered');
Start=min(StackPosition_registered);
End=max(StackPosition_registered+Sizes-1);
Shift=1-Start;
StackRegistered_Translation=zeros(End-Shift,class(IM));
for i = 1:size(StackList,1)
    StackRegistered_Translation(StackPosition_registered(i,1)+Shift(1):StackPosition_registered(i,1)+Shift(1)+Sizes(i,1)-1,...
        StackPosition_registered(i,2)+Shift(2):StackPosition_registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
end



% Affine
for i = 1:size(StackList,1)
    IM = imread(char(StackList(i,1)));
%     Stack(:,:,i) = IM; 
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


for i = 2:size(StackList,1)
    cor3 = corrcoef(double(Stack(:,:,i)),double(Stack(:,:,i-1)));
    correl_before(i-1) = cor3(1,2);
    % Translation
    cor3 = corrcoef(double(StackRegistered_Translation(:,:,i)),double(StackRegistered_Translation(:,:,i-1)));
    correl_Translation(i-1) = cor3(1,2);
    % Affine
    cor3 = corrcoef(double(StackRegistered_Affine(:,:,i)),double(StackRegistered_Affine(:,:,i-1)));
    correl_Affine(i-1) = cor3(1,2);
end

mean(correl_before)
mean(correl_Translation)
mean(correl_Affine)

figure,
plot(correl_before)
hold on
plot(correl_Translation)
plot(correl_Affine)

figure,
plot(sort(correl_before))
hold on
plot(sort(correl_Translation))
plot(sort(correl_Affine))
% for i = 1:size(StackList,1)
%     IM = imread(char(StackList(i,1)));
%     Stack(:,:,i) = IM; 
%     IM = imtranslate(IM,T_Translation.T.b([2,1],i)');
%     StackRegistered(:,:,i) = IM; 
% end
% figure,imshow(max(Stack,[],3),[0 max(Stack(:))]);
% figure,imshow(max(StackRegistered,[],3),[0 max(StackRegistered(:))]);

figure,imshow(max(Stack,[],3),[0 100]);
figure,imshow(max(StackRegistered_Translation,[],3),[0 100]);
figure,imshow(max(StackRegistered_Affine,[],3),[0 100]);




