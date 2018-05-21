% This function finds optimal linear transformations of multiple setts of points (X -> X') 
% X is a cell aray of matched positions in which X{i,j} is [X_i;X_j], 6x(# matches)
% Transform_Type has to be 'Translation', 'Rigid', or 'Affine'
% T contains transformation

function [SaveLocation,StackPositions_Registered,T]=Global_Linear_Transform(StackPositions,FeaturePositions,Transform_Type,DataFolder)

% StackPositions=StackPositions(1:2,:);
% FeaturePositions=FeaturePositions(1:2,1:2);
addpath('../Functions');
parameters;
SaveLocationlatest = [DataFolder,params.G.StackPositions_Registered];

if strcmp(Transform_Type,'Translation')
    delta=0.1; % overal regularization of shifts
    lambda=0; % individual regularization of shifts
    mu=0*10^-6;% regularization of L1 norm
    learning_rate=0.1;
    L=[];
    SaveLocation = [DataFolder,params.GT.StackPositions_Registered];
    SaveLocationCSV = [DataFolder,params.GT.StackPositions_RegisteredCSV];
    TSaveLocation = [DataFolder,'\T_Translation.mat'];
elseif strcmp(Transform_Type,'Rigid')
    alpha=0.1; % overall regularization on Ls
    betta=0; % individual regularizations on Ls
    gamma=1; % enforcing unitary transforms 
    delta=0.1; % overal regularization of shifts
    lambda=0; % individual regularization of shifts
    mu=10;% regularization of L1 norm
    learning_rate=0.1;
    SaveLocation = [DataFolder,params.GR.StackPositions_Registered];
    SaveLocationCSV = [DataFolder,params.GR.StackPositions_RegisteredCSV];
    TSaveLocation = [DataFolder,'\T_Rigid.mat'];
elseif strcmp(Transform_Type,'Affine')
    alpha=0.1; % overall regularization on Ls
    betta=0.1; % individual regularizations on Ls
    delta=0.1; % overal regularization of shifts
    lambda=0; % individual regularization of shifts
    mu=0*10^-6;% regularization of L1 norm
    learning_rate=0.1;
    SaveLocation = [DataFolder,params.GA.StackPositions_Registered];
    SaveLocationCSV = [DataFolder,params.GA.StackPositions_RegisteredCSV];
    TSaveLocation = [DataFolder,'\T_Affine.mat'];
else
    error('Transform_Type has to be Translation, Rigid, or Affine')
    SaveLocation = [DataFolder,params.GN.StackPositions_Registered];
    TSaveLocation = [DataFolder,'\T_nonrigid.mat'];
end

N=size(FeaturePositions,1);
X=FeaturePositions;
for i=1:N
    for j=i+1:N
        if ~isempty(X{i,j})
            X{i,j}(:,1:3)=X{i,j}(:,1:3)+ones(size(X{i,j},1),1)*StackPositions(i,:);
            X{i,j}(:,4:6)=X{i,j}(:,4:6)+ones(size(X{i,j},1),1)*StackPositions(j,:);
        end
    end
end
R_temp=cell2mat(X(:));
disp(mean(sum((R_temp(:,1:3)-R_temp(:,4:6)).^2,2)).^0.5);
R_temp=[R_temp(:,1:3);R_temp(:,4:6)];
Mean=mean(R_temp,1);
Std=std(R_temp,[],1);

for i=1:N
    for j=i+1:N
        if ~isempty(X{i,j})
            X{i,j}=(X{i,j}-ones(size(X{i,j},1),1)*[Mean,Mean])./(ones(size(X{i,j},1),1)*[Std,Std]);
            X{i,j}=X{i,j}';
            X{j,i}=[X{i,j}(4:6,:);X{i,j}(1:3,:)];
        end
    end
end

k=zeros(N,N);
for i=1:N
    for j=1:N
        if ~isempty(X{i,j})
            k(i,j)=size(X{i,j},2);
        end
    end
end
K=k-diag(sum(k));

sum_r=zeros(3*N,N);
Sum_r=zeros(3*N,N);
Sum_rr=zeros(3*N,3*N);
for i=1:N
    temp=zeros(3,3);
    for j=1:N
        if ~isempty(X{i,j})
            sum_r(1+3*(i-1):3*i,j)=sum(X{i,j}(1:3,:),2);
            Sum_rr(1+3*(i-1):3*i,1+3*(j-1):3*j)=X{i,j}(1:3,:)*X{i,j}(4:6,:)';
            temp=temp+X{i,j}(1:3,:)*X{i,j}(1:3,:)';
        end
    end
    Sum_r(1+3*(i-1):3*i,:)=sum_r(1+3*(i-1):3*i,:);
    Sum_r(1+3*(i-1):3*i,i)=sum_r(1+3*(i-1):3*i,i)-sum(sum_r(1+3*(i-1):3*i,:),2);
    Sum_rr(1+3*(i-1):3*i,1+3*(i-1):3*i)=Sum_rr(1+3*(i-1):3*i,1+3*(i-1):3*i)-temp;
end

dEds=zeros(3*N,N);
ss=zeros(3,3*N);
if strcmp(Transform_Type,'Translation')
    I=diag([1,1,1]);
    IN=repmat(I,1,N);
    b=zeros(3,N);
    %s=zeros(3*N,N);
    for i=1:10000
        dEdb=-4/sum(k(:)).*(IN*Sum_r+b*K)+2*delta/N^2.*repmat(sum(b,2),1,N)+2*lambda/N.*b;
        b=b-learning_rate.*dEdb;
    end
elseif strcmp(Transform_Type,'Rigid')
    I=diag([1,1,1]);
    IN=repmat(I,1,N);
    L=IN;
    b=zeros(3,N);
    s=zeros(3*N,N);
    LLL=zeros(size(L));
    for i=1:10000
        for j=1:N
            temp=L(:,1+3*(j-1):3*j);
            LLL(:,1+3*(j-1):3*j)=(temp*temp'-I)*temp;
        end
        ss(1,1:3:end)=0.5.*(diag(Sum_r(1:3:end,:)*(s(1:3:end,:)-s(1:3:end,:)')));
        ss(1,2:3:end)=0.5.*(diag(Sum_r(2:3:end,:)*(s(1:3:end,:)-s(1:3:end,:)')));
        ss(1,3:3:end)=0.5.*(diag(Sum_r(3:3:end,:)*(s(1:3:end,:)-s(1:3:end,:)')));
        
        ss(2,1:3:end)=0.5.*(diag(Sum_r(1:3:end,:)*(s(2:3:end,:)-s(2:3:end,:)')));
        ss(2,2:3:end)=0.5.*(diag(Sum_r(2:3:end,:)*(s(2:3:end,:)-s(2:3:end,:)')));
        ss(2,3:3:end)=0.5.*(diag(Sum_r(3:3:end,:)*(s(2:3:end,:)-s(2:3:end,:)')));
        
        ss(3,1:3:end)=0.5.*(diag(Sum_r(1:3:end,:)*(s(3:3:end,:)-s(3:3:end,:)')));
        ss(3,2:3:end)=0.5.*(diag(Sum_r(2:3:end,:)*(s(3:3:end,:)-s(3:3:end,:)')));
        ss(3,3:3:end)=0.5.*(diag(Sum_r(3:3:end,:)*(s(3:3:end,:)-s(3:3:end,:)')));
        

        dEdL=-4/sum(k(:)).*(L*Sum_rr+b*Sum_r'+ss)+2*alpha/N.*repmat((L*IN'./N-I),1,N)+2*betta/N.*(L-IN)+4*gamma/N.*LLL;
        dEdb=-4/sum(k(:)).*(L*Sum_r+b*K+0.5.*(IN*s-reshape(sum(s,2),3,N)))+2*delta/N^2.*repmat(sum(b,2),1,N)+2*lambda/N.*b;
        
        for p=1:N
            for q=1:N
%                 temp=s(1+3*(p-1):3*p,q);
%                 temp2=zeros(size(temp));
%                 temp2(temp>mu*k(p,q)/sum(k(:)))=temp(temp>mu*k(p,q)/sum(k(:)))-mu*k(p,q)/sum(k(:));
%                 temp2(temp<-mu*k(p,q)/sum(k(:)))=temp(temp<-mu*k(p,q)/sum(k(:)))+mu*k(p,q)/sum(k(:));
                dEds(1+3*(p-1):3*p,q)=2/sum(k(:)).*(L(:,1+3*(p-1):3*p)*sum_r(1+3*(p-1):3*p,q)-L(:,1+3*(q-1):3*q)*sum_r(1+3*(q-1):3*q,p)+k(p,q).*(b(:,p)-b(:,q)+s(1+3*(p-1):3*p,q)))+... 
                    mu.*k(p,q)/sum(k(:))*sign(s(1+3*(p-1):3*p,q));
                    %temp2;
            end
        end
        L=L-learning_rate.*dEdL;
        b=b-learning_rate.*dEdb;
        %s=s-learning_rate.*dEds;
        %s(abs(s)<0.001)=0;
    end    
elseif strcmp(Transform_Type,'Affine')
    I=diag([1,1,1]);
    IN=repmat(I,1,N);
    L=IN;
    b=zeros(3,N);
    %s=zeros(3*N,N);
    for i=1:10000
        dEdL=-4/sum(k(:)).*(L*Sum_rr+b*Sum_r')+2*alpha/N.*repmat((L*IN'./N-I),1,N)+2*betta/N.*(L-IN);
        dEdb=-4/sum(k(:)).*(L*Sum_r+b*K)+2*delta/N^2.*repmat(sum(b,2),1,N)+2*lambda/N.*b;
        L=L-learning_rate.*dEdL;
        b=b-learning_rate.*dEdb;
    end    
end

X_aligned=FeaturePositions;
for i=1:N
    for j=i+1:N
        if k(i,j)>0
            if strcmp(Transform_Type,'Translation')
                X_aligned{i,j}=[X{i,j}(1:3,:)+b(:,i)*ones(1,k(i,j)); X{i,j}(4:6,:)+b(:,j)*ones(1,k(i,j))]';
            elseif strcmp(Transform_Type,'Rigid') || strcmp(Transform_Type,'Affine')
                X_aligned{i,j}=[L(:,1+3*(i-1):3*i)*X{i,j}(1:3,:)+b(:,i)*ones(1,k(i,j));...
                    L(:,1+3*(j-1):3*j)*X{i,j}(4:6,:)+b(:,j)*ones(1,k(i,j))]';
            end
            X_aligned{i,j}=X_aligned{i,j}.*(ones(size(X_aligned{i,j},1),1)*[Std,Std])+ones(size(X_aligned{i,j},1),1)*[Mean,Mean];
        end
    end
end

for i=1:N
    if strcmp(Transform_Type,'Translation')
        b(:,i)=b(:,i).*Std';
    elseif strcmp(Transform_Type,'Rigid') || strcmp(Transform_Type,'Affine')
        b(:,i)=b(:,i).*Std'+Mean'-L(:,1+3*(i-1):3*i)*Mean';
    end
end

% if strcmp(Transform_Type,'Translation') 
%     StackPositions_Registered =  (StackPositions(:,:)) + b(:,:)';
% elseif strcmp(Transform_Type,'Rigid') || strcmp(Transform_Type,'Affine')
%     StackPositions_Registered(i,j)=[L(:,1+3*(i-1):3*i)*StackPositions(i,j)(1:3,:)+b(:,i)*ones(1,k(i,j));...
%         L(:,1+3*(j-1):3*j)*StackPositions(i,j)(4:6,:)+b(:,j)*ones(1,k(i,j))]';  
% end

StackPositions_Registered =  (StackPositions(:,:)) + b(:,:)';

T.transform = Transform_Type;
T.b = b;
T.L = L;
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
save(SaveLocationlatest,'StackPositions_Registered');
save(SaveLocation,'StackPositions_Registered');
csvwrite(SaveLocationCSV,StackPositions_Registered);
csvwrite([DataFolder,params.G.StackPositions_RegisteredCSV],StackPositions_Registered);
disp(['Registered Positions Data saved as: ',SaveLocationlatest]);

R_temp=cell2mat(X_aligned(:));
disp(mean(sum((R_temp(:,1:3)-R_temp(:,4:6)).^2,2)).^0.5);
%sum(sum(abs(s(1:3:end,:))+abs(s(2:3:end,:))+abs(s(3:3:end,:)).*k))


