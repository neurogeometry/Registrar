% This function finds optimal translations of multiple setts of points (X -> X') 
% X is a cell aray of matched positions in which X{i,j} is [X_i,X_j], (# matches)x6

function [SaveLocation,T]=Global_Optimal_Translation(StackPositions_pixels,X,DataFolder)
addpath('../Functions');
parameters;
SaveLocation = [DataFolder,params.GT.StackPositions_Registered];
TSaveLocation = [DataFolder,'\T_Translation.mat'];
% regularization of the overall shift 
eta=params.GT.eta;
% regularization of individual shifts
lambda=params.GT.lamda;

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

T.transform = 'Translation';
T.b = b;
T.L = [];
T.N_L = [];
T.P = [];
T.Min = [];
T.Max = [];
T.XYZlmn = [];

% StackPositions_Registered=cell(N,N);
% for i=1:N
%     for j=i+1:N
%         if k(i,j)>0
%             StackPositions_Registered{i,j}=StackPositions_pixels(i,j)+ones(k(i,j),1)*[b(:,i);b(:,j)]';
%         end
%     end
% end
save(TSaveLocation,'T');
save(SaveLocation,'StackPositions_Registered');
csvwrite([DataFolder,params.GT.StackPositions_RegisteredCSV],StackPositions_Registered);
disp(['Registered Positions Data saved as: ',SaveLocation]);
