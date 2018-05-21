function [SaveLocation,StackPositions_Registered]=GlobalOptimization(StackPositions_pixels,Matched,DataFolder)
addpath('../Functions');
parameters;
SaveLocation = [DataFolder,params.GO.StackPositions_Registered];
GlobalMatchd=Matched;
for i = 1:size(Matched,1)
    for j=1:size(Matched,2)
        if i<j && ~isempty(Matched{i,j})
            GlobalMatchd{i,j}=[Matched{i,j}(:,1:3)+StackPositions_pixels(i,:)-1,Matched{i,j}(:,4:6)+StackPositions_pixels(j,:)-1];
        end
        if i>j && ~isempty(Matched{j,i})
            GlobalMatchd{i,j}=[Matched{j,i}(:,4:6)+StackPositions_pixels(i,:)-1,Matched{j,i}(:,1:3)+StackPositions_pixels(j,:)-1];
        end
    end
end


K=zeros(size(GlobalMatchd));
R=zeros(size(GlobalMatchd,1),3);
for i = 1:size(GlobalMatchd,1)
    Rtemp=zeros(1,3);
    for j=1:size(GlobalMatchd,2)
        M=GlobalMatchd{i,j};
        if ~isempty(M)
            K(i,j)=size(M,1);
            Rtemp=Rtemp+sum(M(:,1:3)-M(:,4:6),1);
        end
    end
    R(i,:)=Rtemp;
end
K=K-diag(sum(K)+params.GO.eps)-params.GO.lamda/2;
b=K\R;

% StackPositions_Registered_NEW =  (StackPositions_pixels(:,:)) + b(:,:).*resolution.*1000
StackPositions_Registered =  (StackPositions_pixels(:,:)) + b(:,:);

% sourceID = 21;
% targetID = 20;
% StackPositions_Registered_NEW(sourceID,:) - StackPositions_Registered_NEW(targetID,:)
% StackPositions_Registered_NEW1(sourceID,:) - StackPositions_Registered_NEW1(targetID,:)
% StackPositions_pixels(sourceID,:)-StackPositions_pixels(targetID,:)

save(SaveLocation,'StackPositions_Registered');
csvwrite([DataFolder,params.GO.StackPositions_RegisteredCSV],StackPositions_Registered);
disp(['Registered Positions Data saved as: ',SaveLocation]);
end