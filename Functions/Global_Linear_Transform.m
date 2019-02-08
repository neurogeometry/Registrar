% This function finds optimal linear transformations of multiple setts of points (X -> X')
% X is a cell aray of matched positions in which X{i,j} is [X_i;X_j], 6x(# matches)
% Transform_Type has to be 'Translation', 'Rigid', or 'Affine'
% Backward stepwise selection is used to eliminate outliers in matches
% T contains transformation

function [TSaveLocation,T]=Global_Linear_Transform(StackPositions,FeaturePositions,Transform_Type,DataFolder)
addpath('../Functions');
parameters;
SaveLocationlatest = [DataFolder,params.G.StackPositions_Registered];

if strcmp(Transform_Type,'Translation')
    delta=0.1; % overal regularization of shifts
    lambda=0; % individual regularization of shifts
    learning_rate=0.1;
    Max_Iterations=10000;
    ObjFunc_tol=10^-8;
    Dist_tol=2; % pixels
    Conn_tol=2;
    n_dist=10;
    
    L=[];
    TSaveLocation = [DataFolder,'/T_Translation.mat'];
elseif strcmp(Transform_Type,'Rigid')
    alpha=0.1; % overall regularization on Ls
    betta=0.1; % individual regularizations on Ls
    gamma=5; % enforcing unitary transforms
    delta=0.1; % overal regularization of shifts
    lambda=0; % individual regularization of shifts
    learning_rate=0.1;
    Max_Iterations=10000;
    ObjFunc_tol=10^-8;
    Dist_tol=2; % pixels
    Conn_tol=3;
    n_dist=10;
    TSaveLocation = [DataFolder,'/T_Rigid.mat'];
elseif strcmp(Transform_Type,'Affine')
    alpha=0.1; % overall regularization on Ls
    betta=0.1; % individual regularizations on Ls
    delta=0.1; % overal regularization of shifts
    lambda=0; % individual regularization of shifts
    learning_rate=0.1;
    Max_Iterations=10000;
    ObjFunc_tol=10^-8;
    Dist_tol=2; % pixels
    Conn_tol=3;
    n_dist=10;
%     SaveLocation = [DataFolder,params.GA.StackPositions_Registered];
    TSaveLocation = [DataFolder,'/T_Affine.mat'];
else
    error('Transform_Type has to be Translation, Rigid, or Affine')
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
Std=max(std(R_temp,[],1)).*[1,1,1];

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

I=diag([1,1,1]);
IN=repmat(I,1,N);
outlier=true;
while outlier
    K=k-diag(sum(k));
    
    sum_r=zeros(3*N,N);
    Sum_r=zeros(3*N,N);
    Sum_rr=zeros(3*N,3*N);
    for i=1:N
        temp=zeros(3,3);
        for j=1:N
            if k(i,j)>0
                sum_r(1+3*(i-1):3*i,j)=sum(X{i,j}(1:3,:),2);
                Sum_rr(1+3*(i-1):3*i,1+3*(j-1):3*j)=X{i,j}(1:3,:)*X{i,j}(4:6,:)';
                temp=temp+X{i,j}(1:3,:)*X{i,j}(1:3,:)';
            end
        end
        Sum_r(1+3*(i-1):3*i,:)=sum_r(1+3*(i-1):3*i,:);
        Sum_r(1+3*(i-1):3*i,i)=sum_r(1+3*(i-1):3*i,i)-sum(sum_r(1+3*(i-1):3*i,:),2);
        Sum_rr(1+3*(i-1):3*i,1+3*(i-1):3*i)=Sum_rr(1+3*(i-1):3*i,1+3*(i-1):3*i)-temp;
    end
    
    i=0;
    RMS_dist=zeros(1,ceil(Max_Iterations/n_dist));
    flag=true;
    if strcmp(Transform_Type,'Translation')
        b=zeros(3,N);
        while i<=Max_Iterations && flag
            i=i+1;
            dEdb=-4/sum(k(:)).*(IN*Sum_r+b*K)+2*delta/N^2.*repmat(sum(b,2),1,N)+2*lambda/N.*b;
            b=b-learning_rate.*dEdb;
            
            if mod(i,n_dist)==0
                RMS_temp=0;
                for ii=1:N
                    for jj=ii+1:N
                        if k(ii,jj)>0
                            RMS_temp=RMS_temp+sum(sum((X{ii,jj}(1:3,:)+b(:,ii)*ones(1,k(ii,jj))-X{ii,jj}(4:6,:)-b(:,jj)*ones(1,k(ii,jj))).^2));
                        end
                    end
                end
                RMS_dist(i/n_dist)=(RMS_temp/(sum(k(:))/2))^0.5;
                if i/n_dist>1
                    flag=abs(RMS_dist(i/n_dist)-RMS_dist(i/n_dist-1))/(RMS_dist(i/n_dist)+RMS_dist(i/n_dist-1))>ObjFunc_tol;
                end
            end
        end
    elseif strcmp(Transform_Type,'Rigid')
        L=IN;
        b=zeros(3,N);
        LLL=zeros(size(L));
        while i<=Max_Iterations && flag
            i=i+1;
            for j=1:N
                temp=L(:,1+3*(j-1):3*j);
                LLL(:,1+3*(j-1):3*j)=(temp*temp'-I)*temp;
            end
            dEdL=-4/sum(k(:)).*(L*Sum_rr+b*Sum_r')+2*alpha/N.*repmat((L*IN'./N-I),1,N)+2*betta/N.*(L-IN)+4*gamma/N.*LLL;
            dEdb=-4/sum(k(:)).*(L*Sum_r+b*K)+2*delta/N^2.*repmat(sum(b,2),1,N)+2*lambda/N.*b;
            L=L-learning_rate.*dEdL;
            b=b-learning_rate.*dEdb;
            
            if mod(i,n_dist)==0
                RMS_temp=0;
                for ii=1:N
                    for jj=ii+1:N
                        if k(ii,jj)>0
                            RMS_temp=RMS_temp+sum(sum((L(:,1+3*(ii-1):3*ii)*X{ii,jj}(1:3,:)+b(:,ii)*ones(1,k(ii,jj))-L(:,1+3*(jj-1):3*jj)*X{ii,jj}(4:6,:)-b(:,jj)*ones(1,k(ii,jj))).^2));
                        end
                    end
                end
                RMS_dist(i/n_dist)=(RMS_temp/(sum(k(:))/2))^0.5;
                if i/n_dist>1
                    flag=abs(RMS_dist(i/n_dist)-RMS_dist(i/n_dist-1))/(RMS_dist(i/n_dist)+RMS_dist(i/n_dist-1))>ObjFunc_tol;
                end
            end
        end
    elseif strcmp(Transform_Type,'Affine')
        L=IN;
        b=zeros(3,N);
        while i<=Max_Iterations && flag
            i=i+1;
            dEdL=-4/sum(k(:)).*(L*Sum_rr+b*Sum_r')+2*alpha/N.*repmat((L*IN'./N-I),1,N)+2*betta/N.*(L-IN);
            dEdb=-4/sum(k(:)).*(L*Sum_r+b*K)+2*delta/N^2.*repmat(sum(b,2),1,N)+2*lambda/N.*b;
            L=L-learning_rate.*dEdL;
            b=b-learning_rate.*dEdb;
            
            if mod(i,n_dist)==0
                RMS_temp=0;
                for ii=1:N
                    for jj=ii+1:N
                        if k(ii,jj)>0
                            RMS_temp=RMS_temp+sum(sum((L(:,1+3*(ii-1):3*ii)*X{ii,jj}(1:3,:)+b(:,ii)*ones(1,k(ii,jj))-L(:,1+3*(jj-1):3*jj)*X{ii,jj}(4:6,:)-b(:,jj)*ones(1,k(ii,jj))).^2));
                        end
                    end
                end
                RMS_dist(i/n_dist)=(RMS_temp/(sum(k(:))/2))^0.5;
                if i/n_dist>1
                    flag=abs(RMS_dist(i/n_dist)-RMS_dist(i/n_dist-1))/(RMS_dist(i/n_dist)+RMS_dist(i/n_dist-1))>ObjFunc_tol;
                end
            end
        end
    end
    
    RMS_final=zeros(N,N);
    for ii=1:N
        for jj=ii+1:N
            if k(ii,jj)>0
                if strcmp(Transform_Type,'Translation')
                    X_aligned_ii=X{ii,jj}(1:3,:)+b(:,ii)*ones(1,k(ii,jj));
                    X_aligned_jj=X{ii,jj}(4:6,:)+b(:,jj)*ones(1,k(ii,jj));
                elseif strcmp(Transform_Type,'Rigid') || strcmp(Transform_Type,'Affine')
                    X_aligned_ii=L(:,1+3*(ii-1):3*ii)*X{ii,jj}(1:3,:)+b(:,ii)*ones(1,k(ii,jj));
                    X_aligned_jj=L(:,1+3*(jj-1):3*jj)*X{ii,jj}(4:6,:)+b(:,jj)*ones(1,k(ii,jj));
                end
                X_aligned_ii=X_aligned_ii.*(Std'*ones(1,size(X_aligned_ii,2)))+Mean'*ones(1,size(X_aligned_ii,2));
                X_aligned_jj=X_aligned_jj.*(Std'*ones(1,size(X_aligned_jj,2)))+Mean'*ones(1,size(X_aligned_jj,2));
                RMS_final(ii,jj)=(mean(sum((X_aligned_ii-X_aligned_jj).^2,1)))^0.5;
            end
        end
    end
    
    M=max(RMS_final(:));
    [iii,jjj]=find(RMS_final==M);
    if M>=Dist_tol && k(iii,jjj)<=Conn_tol
        k(iii,jjj)=0; k(jjj,iii)=0;
        disp('Connection Removed')
    else
        outlier=false;
    end
end


X_aligned=cell(N,N);
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
R_temp=cell2mat(X_aligned(:));
disp(mean(sum((R_temp(:,1:3)-R_temp(:,4:6)).^2,2)).^0.5);


for i=1:N
    if strcmp(Transform_Type,'Translation')
        b(:,i)=b(:,i).*Std';
    elseif strcmp(Transform_Type,'Rigid') || strcmp(Transform_Type,'Affine')
        b(:,i)=b(:,i).*Std'+Mean'-L(:,1+3*(i-1):3*i)*Mean';
    end
end

T.transform=Transform_Type;
T.b=b;
T.L=L;
T.N_L=[];
T.P=[];
T.Min=[];
T.Max=[];
T.XYZlmn=[];

save(TSaveLocation,'T');
disp(['Registered Positions Data saved as: ',TSaveLocation]);