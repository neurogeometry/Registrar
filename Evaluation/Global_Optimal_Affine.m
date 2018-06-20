% This function finds optimal affine transformations of multiple setts of points (X -> X') 
% X is a cell aray of matched positions in which X{i,j} is [X_i;X_j], 6x(# matches)

function [X_aligned,L,b]=Global_Optimal_Affine(StackPositions,FeaturePositions)


% StackPositions=StackPositions(1:2,:);
% FeaturePositions=FeaturePositions(1:2,1:2);

% regularization of the overall shift 
eta=1;
% regularization of individual shifts
lambda=10^-6;
% regularization of the overall deformation 
nu=1;%10^6;
% regularization of individual deformations 
mu=10^6;

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
% disp(mean(sum((R_temp(:,1:3)-R_temp(:,4:6)).^2,2)).^0.5);
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
K=k-diag(sum(k))-(lambda/2/N*sum(k(:))).*diag(ones(1,N))-(eta/2/N^2).*sum(k(:));

r=zeros(3*N,N);
R=zeros(3*N,3*N);
for i=1:N
    temp=zeros(3,3);
    for j=1:N
        if ~isempty(X{i,j})
            r(1+3*(i-1):3*i,j)=sum(X{i,j}(1:3,:),2);
            R(1+3*(i-1):3*i,1+3*(j-1):3*j)=X{i,j}(1:3,:)*X{i,j}(4:6,:)';
            temp=temp+X{i,j}(1:3,:)*X{i,j}(1:3,:)';
        end
    end
    r(1+3*(i-1):3*i,i)=r(1+3*(i-1):3*i,i)-sum(r(1+3*(i-1):3*i,:),2);
    R(1+3*(i-1):3*i,1+3*(i-1):3*i)=R(1+3*(i-1):3*i,1+3*(i-1):3*i)-temp-(mu/2/N*sum(k(:))).*diag([1,1,1]);
    R(1+3*(i-1):3*i,:)=R(1+3*(i-1):3*i,:)-(nu/2/N^2*sum(k(:))).*repmat(diag([1,1,1]),1,N);
end
I=((mu+nu)/2/N*sum(k(:))).*repmat(diag([1,1,1]),1,N);
L=I/((r/K)*r'-R);
b=-L*r/K;

X_aligned=cell(N,N);
for i=1:N
    for j=i+1:N
        if k(i,j)>0
            X_aligned{i,j}=[L(:,1+3*(i-1):3*i)*X{i,j}(1:3,:)+b(:,i)*ones(1,k(i,j));...
                L(:,1+3*(j-1):3*j)*X{i,j}(4:6,:)+b(:,j)*ones(1,k(i,j))]';
            X_aligned{i,j}=X_aligned{i,j}.*(ones(size(X_aligned{i,j},1),1)*[Std,Std])+ones(size(X_aligned{i,j},1),1)*[Mean,Mean];
        end
    end
end

for i=1:N
    b(:,i)=b(:,i).*Std'+Mean'-L(:,1+3*(i-1):3*i)*Mean';
end

R_temp=cell2mat(X_aligned(:));
% disp(mean(sum((R_temp(:,1:3)-R_temp(:,4:6)).^2,2)).^0.5);


