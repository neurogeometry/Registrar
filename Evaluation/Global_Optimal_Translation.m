function [SaveLocation,StackPositions_Registered]=Global_Optimal_Translation(StackPositions_pixels,X,DataFolder)
addpath('../Functions');
parameters;
SaveLocation = [DataFolder,params.GO.StackPositions_Registered];

% regularization of the overall shift 
eta=1;
lambda=10^-6;

N=size(X,1);

for i=1:N
    for j=i+1:N
        if ~isempty(X{i,j})
            X{i,j}(:,1:3)=X{i,j}(:,1:3)+ones(size(X{i,j},1),1)*StackPositions_pixels(i,:)-1;
            X{i,j}(:,4:6)=X{i,j}(:,4:6)+ones(size(X{i,j},1),1)*StackPositions_pixels(j,:)-1;
            X{i,j}=X{i,j}';
            X{j,i}=[X{i,j}(4:6,:);X{i,j}(1:3,:)];
        end
    end
end

del_r=zeros(3,N);
k=zeros(N,N);

for i=1:N
    temp=zeros(3,N);
    for j=1:N
        if ~isempty(X{i,j})
            k(i,j)=size(X{i,j},2);
            temp(:,j)=sum(X{i,j}(1:3,:)-X{i,j}(4:6,:),2);
        end
    end
    del_r(:,i)=sum(temp,2);
end
K=k-diag(sum(k))-eta/2/N^2*sum(k(:))-lambda/2/N*sum(k(:))*diag(ones(1,N));
b=del_r/K;

StackPositions_Registered =  (StackPositions_pixels(:,:)) + b(:,:)';

save(SaveLocation,'StackPositions_Registered');
csvwrite([DataFolder,params.GO.StackPositions_RegisteredCSV],StackPositions_Registered);
disp(['Registered Positions Data saved as: ',SaveLocation]);
