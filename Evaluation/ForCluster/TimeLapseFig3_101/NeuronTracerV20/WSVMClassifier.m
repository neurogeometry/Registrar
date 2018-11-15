% This function performs classification based on a weighted linear kernel svm algorithm.
% In each group, cost of the correct merger has to be lower than all the rest.
% SVM is modified to work with 1 class (+1) and no bias.
% Weights reflect user confidence in the merged examples.
% This version is designed to work with ClusterStrW and uses previous TrainingHistory as an input

function ClustersStr = WSVMClassifier(ClustersStr,TrainingHistory,Parameters)

T=Parameters.Classification.ConfidenceScale; % temperature used to define the confidence measure 
C=Parameters.Classification.SVMMultiplierC; % Lagrange multiplier. Inversely correlated with the margin. 
stdX=[80;12;60;9;3;3;0.5;0.2;12];

disp('Weighted SVM classification started.')

if isempty(TrainingHistory)
    X_0={ClustersStr.cost_components};
    Xp_0={ClustersStr.best_merger};
    alpha_0={ClustersStr.alpha};
    q_0=[ClustersStr.weights]; % user defined weights of the best mergers
else
    X_0=[{ClustersStr.cost_components},{TrainingHistory.cost_components}];
    Xp_0=[{ClustersStr.best_merger},{TrainingHistory.best_merger}];
    alpha_0=[{ClustersStr.alpha},{TrainingHistory.alpha}];
    q_0=[[ClustersStr.weights],[TrainingHistory.weights]]; % user defined weights of the best mergers
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extra dimentions and subtruction
X=X_0;
q=cell(size(Xp_0));
for i=1:length(X_0)  
    X{i}=X{i}./(stdX*ones(1,size(X{i},2))); % normalization    
    q{i}=q_0(i).*ones(1,length(Xp_0{i}));
    temp=find(Xp_0{i}==1);
    if ~isempty(temp) && q_0(i)>0 % training is only done on weighted examples
        Xp_0{i}(temp)=NaN;
        X{i}=X{i}-X{i}(:,temp)*ones(1,length(Xp_0{i}));
    else
        Xp_0{i}=NaN.*Xp_0{i};
    end
end

X=cell2mat(X); 
Xp_0=cell2mat(Xp_0);
alpha=cell2mat(alpha_0);
q=cell2mat(q);

ind=isnan(Xp_0);
q(ind)=[];
alpha(ind)=[];
X(:,ind)=[]; 
X=sparse(X);
alpha=sparse(alpha);
num_cost_terms=size(X,1); % N
num_examples=size(X,2); % m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximization of Ld
if ~isempty(X)
    f=sparse(-ones(length(alpha),1));
    warning('off','optim:quadprog:HessianNotSym')
    
    H=X'*X;
    options=optimset('Algorithm','interior-point-convex','Display','iter-detailed');
    alpha=quadprog(H./num_examples^2./num_cost_terms,f/num_examples,[],[],[],[],zeros(length(alpha),1),C.*q',[],options);
%     alpha=quadprog(H./num_examples^2./num_cost_terms,f/num_examples,[],[],[],[],zeros(length(alpha),1),C.*q',alpha,options);
    
%     options=optimset('Algorithm','trust-region-reflective','HessMult',@qpbox4mult,'TolFun',10^-8,'TolX',10^-8,'Display','iter-detailed'); % 
%     alpha=quadprog(X./(num_examples^2*num_cost_terms)^0.5,f/num_examples,[],[],[],[],zeros(length(alpha),1),C.*q',[],options);
        
%    alpha(alpha<10^-8)=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha=alpha';
else
    alpha=[];
end

%temp=(alpha(alpha>0)*Xprime(alpha>0,:))*X;
%i_F=nnz(temp>=0);
%Ld=sum(alpha)-temp*alpha'/2; % dual Lagrangian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% support vectors 
svs=X(:,alpha>0);
if ~isempty(svs)
    % weights for normalized and original cost components
    w_norm=alpha(alpha>0)*svs'./num_examples;
    w_orig=w_norm./stdX';
        
    display(['M_norm = ' num2str(1/sum(w_norm.^2)^0.5)])
    display(['M_orig = ' num2str(1/sum(w_orig.^2)^0.5)])
    
    %w_orig=w_orig./sum(w_orig.^2)^0.5;    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update classifications, Xp
    Xp=cell(size(X_0));
    alpha_temp=zeros(size(ind));
    alpha_temp(~ind)=alpha;
    count=0;
    Cluster_errors=0;
    Cluster_number=0;
    Scenario_errors=0;
    Scenario_number=0;
    for i=1:size(ClustersStr,2)
        ClustersStr(i).w=w_orig;
                
        Scenario_costs=w_orig*X_0{i};
        [min_Cost,Best_detected_scanario]=min(Scenario_costs);
        Best_labeled_scanario=find(ClustersStr(i).best_merger==1);
        
        Xp{i}=-ones(1,size(X_0{i},2));
        Xp{i}(Best_detected_scanario)=1;
        ClustersStr(i).confidence=1/sum(exp((min_Cost-Scenario_costs)./T));
        
        if ~isempty(Best_labeled_scanario) && q_0(i)>0 % errors are only calculated for labeled and weighted examples
            Cluster_number=Cluster_number+1;
            Scenario_number=Scenario_number+length(Scenario_costs);
                      
            if nnz(ClustersStr(i).best_merger~=Xp{i})
                Cluster_errors=Cluster_errors+1;
                Scenario_errors=Scenario_errors+sum(Scenario_costs<Scenario_costs(Best_labeled_scanario));
            end
        end
        ClustersStr(i).best_merger=Xp{i};
        ClustersStr(i).alpha=alpha_temp(count+1:count+size(X_0{i},2));
        count=count+size(X_0{i},2);
    end
    display(['Erroneously merged weighted clusters: ' num2str(Cluster_errors), ' out of ', num2str(Cluster_number)])
    display(['Total number of SVM classification errors: ' num2str(Scenario_errors), ' out of ', num2str(Scenario_number)])
else
    display('There are no valid classification examples.')
end

disp('Weighted SVM classification is complete.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function W = qpbox4mult(A,Y)
% W = A'*(A*Y); %+10^-8.*diag(sparse(ones(1,size(A,2))))*Y;
    

