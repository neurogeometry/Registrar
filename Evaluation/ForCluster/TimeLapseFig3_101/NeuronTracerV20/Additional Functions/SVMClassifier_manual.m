% This function performs classification based on modified linear kernel-based svm algorithm.
% In each group, cost of the correct merger has to be lower than all the rest.
% SVM is modified to work with 1 class (+1) and no bias.
% This version works with X_Xp.mat files

function [Xp,alpha]=SVMClassifier_manual(X_0,Xp_0,alpha_0)

C=0.1; % Lagrange multiplayer. Inversely correlated with the margin. 
Nsteps=1000000; % Number of gradient steps
betta0=0.000003; % gradiat step size

format short g

num_cost_terms=size(X_0{1},1);
[ii,jj]=find(triu(ones(num_cost_terms,num_cost_terms)));

X=X_0;
for i=1:length(X_0)
    Xp_0{i}(sum(X{i})==0)=NaN;
    Xp_0{i}(sum(X{i})==Inf)=NaN;
    Xp_0{i}(isnan(sum(X{i})))=NaN;
    temp=find(Xp_0{i}==1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extra dimentions
    X{i}=[X{i};zeros(length(ii),size(X{i},2))];
    for k=1:length(ii)
        X{i}(num_cost_terms+k,:)=X{i}(ii(k),:).*X{i}(jj(k),:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if length(temp)==1
        Xp_0{i}(temp)=NaN;
        X{i}=X{i}-X{i}(:,temp)*ones(1,length(Xp_0{i}));
    else
        Xp_0{i}=NaN;
        X{i}=nan(size(X{1},1),1);
    end
end

X=cell2mat(X);
Xp_0=cell2mat(Xp_0);
ind=isnan(Xp_0);
X(:,ind)=[]; 
[N,m]=size(X);
stdX=std(X,[],2);
X=X./(stdX*ones(1,m));
Xprime=X';

% kernels must be zero at x=0
%K=X'*X;

if isempty(alpha_0)
    alpha_0=C/2.*ones(1,m);
    %alpha_0=ones(1,m);
end

alpha=alpha_0;
temp=(alpha(alpha>0)*Xprime(alpha>0,:))*X;
Ld=sum(alpha)-temp*alpha'/2; % dual Lagrangian
delta_alpha=1;
Ld_old=-inf;
iteration=1;
while iteration<=Nsteps && (Ld-Ld_old)/(sum(delta_alpha.^2))^0.5>10^-10  
    Ld_old=Ld;
    if iteration>100000
        betta=betta0*10;
    elseif iteration>10000
        betta=betta0*5;
    else
        betta=betta0;
    end
        
    delta_alpha=betta.*(1-temp); % gradient descent
    alpha=alpha+delta_alpha;
    alpha(alpha<0)=0;
    alpha(alpha>C)=C;
    
    %temp=alpha*K;
    %temp=alpha(alpha>0)*K(alpha>0,:);
    temp=(alpha(alpha>0)*Xprime(alpha>0,:))*X;
    i_F=nnz(temp>=0);
    Ld=sum(alpha)-temp*alpha'/2;
    
    if mod(iteration,1000)==0
        [iteration,Ld,(Ld-Ld_old)/(sum(delta_alpha.^2))^0.5,max(abs(delta_alpha)),i_F/m]
    end
    iteration=iteration+1;  
end

svs=X(:,alpha>0);
w=alpha(alpha>0)*svs';
display(['M = ' num2str(2/sum(w.^2)^0.5)])

w=w./stdX';
w=w./sum(w.^2)^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update classification, Xp
Xp=cell(size(X_0));
for i=1:length(X_0)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extra dimentions
    X_0{i}=[X_0{i};zeros(length(ii),size(X_0{i},2))];
    for k=1:length(ii)
        X_0{i}(num_cost_terms+k,:)=X_0{i}(ii(k),:).*X_0{i}(jj(k),:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% figure
% plot([1:iteration],(Ld(1:iteration)))
% figure
% hist(alpha,[0:C])

