% This function performs online classification based on the sign-constrained 
% perceptron learning rule. The classifier attempts to insures that the total 
% costs of the correct merger scenarios within individual end-point clusters are the lowest. 
% Cell array Cost_Components contains cost components for the best (typically top 20) 
% merger scenarios for each end-point cluster. User_Input is an array of two numbers: 
% cluster index and user corrected merger index within that cluster. 

% NOT FINISHED
function [New_w,New_Mergers]=OnlinePerceptronClassifier(Cost_Components,Current_w,User_Input)

kappa=0.001; % weight scaling factor
Nsteps=500000; % Number of training steps

if  isempty(Current_w) && isempty(User_Input)
    New_w=[];
    New_Mergers=[]; % no mergers
end

if  isempty(Current_w) && ~isempty(User_Input)
    
end

for i=1:length(Xp)
    Xp{i}(sum(X{i})==0)=NaN;
    Xp{i}(sum(X{i})==Inf)=NaN;
    Xp{i}(isnan(sum(X{i})))=NaN;
    temp=find(Xp{i}==1);
    if length(temp)==1
        Xp{i}(temp)=NaN;
        X{i}=X{i}-X{i}(:,temp)*ones(1,length(Xp{i}));
    else
        Xp{i}=NaN;
        X{i}=nan(size(X{1},1),1);
    end
end

X=cell2mat(X);
Xp=cell2mat(Xp);
ind=isnan(Xp);
X(:,ind)=[]; 
[N,m]=size(X);
Xp=ones(1,m);
% normalization of X
X=X./(mean(X,2)*ones(1,m));

g=[1,-1,-1,1,1,1,1,1];
if ~isempty(Current_w)
    w=Current_w;
else
    w=g;
end

iteration=0;
delta_w=((ones(N,1)*Xp).*X)'./N^0.5;
out=(kappa/N^0.5-Xp.*(w*X)./N)>0;
learned_associations=m-sum(out);
learned_associations_max=learned_associations;
New_w=w;
while learned_associations<m && iteration<=Nsteps   
    % update of associations one at a time
    %mu0=find(out,1,'first');
    mu0=find(out); mu0=mu0(randi(length(mu0)));
    iteration=iteration+1;
    w=w+delta_w(mu0,:);
    w(w.*g<0)=0;
    out=(kappa/N^0.5-Xp.*(w*X)./N)>0;
    learned_associations=m-sum(out);
    if learned_associations>learned_associations_max
        New_w=w;
        learned_associations_max=learned_associations;
    end
    [iteration,learned_associations,learned_associations_max]
end

if learned_associations_max~=m
    ['There is no feasible solution: ', num2str(learned_associations_max/m*100), '% success rate']
else
    ['Training is complete: ', num2str(learned_associations_max/m*100), '% success rate']
end

%New_w=New_w.*mean(X,2);
New_w=New_w./sum(New_w.^2)^0.5
        
