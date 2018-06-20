% This function performs classification based on perceptron learning rule.
% Costs of the correct merger scenarios in each cluster have to be the lowest. 
% This version works with weights which reflect user-assigned importance of branch merging examples.
% This version is designed to work with ClusterStrW and uses previous TrainingHistory as an input.

function ClustersStr=WPerceptronClassifier(ClustersStr,TrainingHistory,Parameters,varargin)

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

Nsteps=100000; % Number of training steps
kappa=Parameters.Classification.PerceptronMargin; % weight scaling factor
T=Parameters.Classification.ConfidenceScale; % temperature used to define the confidence measure 
stdX=[80;12;60;9;3;3;0.5;0.2;12];

disp('Weighted perceptron classification started.')
disp(['     Iteration    ',   'Learned   ',  'All scenarios'])

if isempty(TrainingHistory)
    X_0={ClustersStr.cost_components};
    Xp_0={ClustersStr.best_merger};
    q_0=[ClustersStr.weights]; % user defined weights of the best mergers
else
    X_0=[{ClustersStr.cost_components},{TrainingHistory.cost_components}];
    Xp_0=[{ClustersStr.best_merger},{TrainingHistory.best_merger}];
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
    if ~isempty(temp) && q_0(i)>0 % training is only done on labeled and weighted examples
        Xp_0{i}(temp)=NaN;
        X{i}=X{i}-X{i}(:,temp)*ones(1,length(Xp_0{i}));
    else
        Xp_0{i}=NaN.*Xp_0{i};
    end
end

X=cell2mat(X);
Xp_0=cell2mat(Xp_0);
q=cell2mat(q);

ind=isnan(Xp_0);
q(ind)=[];
X(:,ind)=[]; 
num_cost_terms=size(X,1); % N
num_examples=size(X,2); % m
Xp_temp=ones(1,num_examples);

% g=[1,-1,-1,1,1,1,1,1];
% w=g;

if isfield(TrainingHistory,'w') && ~isempty(TrainingHistory)
    w=TrainingHistory(1).w;
else
    w=zeros(1,num_cost_terms);
end

if nnz(q)>0
    iteration=0;
    delta_w=((ones(num_cost_terms,1)*Xp_temp).*X)'./num_cost_terms^0.5;
    out=(kappa/num_cost_terms^0.5-Xp_temp.*(w*X)./num_cost_terms)>0;
    learned_associations=num_examples-sum(out);
    learned_associations_max=learned_associations;
    w_max=w;
    while learned_associations<num_examples && iteration<=Nsteps
        if myClass.getFlag()==1
               return;
        end
         
        % update of associations one at a time
        %mu0=find(out,1,'first');
        %mu0=find(out); mu0=mu0(randi(length(mu0)));
        mu0=find(out); 
        q0=q(out); 
        q0=q0./max(q0); 
        mu0=mu0(q0>rand(1,length(q0))); 
        mu0=mu0(randi(length(mu0)));
        
        iteration=iteration+1;
        w=w+delta_w(mu0,:);
        %w(w.*g<0)=0;
        out=(kappa/num_cost_terms^0.5-Xp_temp.*(w*X)./num_cost_terms)>0;
        learned_associations=num_examples-sum(out);
        if learned_associations>learned_associations_max
            w_max=w;
            learned_associations_max=learned_associations;
        end
        if mod(iteration,1000)==0
            disp([iteration,learned_associations_max,num_examples])
        end
    end
    if iteration<=Nsteps
        disp([iteration,learned_associations_max,num_examples])
    end
    
    w_orig=w_max./stdX';
    %w_orig=w_orig./sum(w_orig.^2)^0.5;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update classifications, Xp
    Xp=cell(size(X_0));
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
        count=count+size(X_0{i},2);
    end
    display(['Erroneously merged weighted clusters: ' num2str(Cluster_errors), ' out of ', num2str(Cluster_number)])
    display(['Total number of perceptron classification errors: ' num2str(Scenario_errors), ' out of ', num2str(Scenario_number)])
else
    display('There are no valid classification examples.')
end

disp('Weighted perceptron classification is complete.')