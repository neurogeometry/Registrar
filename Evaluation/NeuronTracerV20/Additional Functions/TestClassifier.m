function TestClassifier(X,Xp,w_perceptron,w_svm)

num_cost_terms=size(X{1},1);
[ii,jj]=find(triu(ones(num_cost_terms,num_cost_terms)));

for i=1:length(Xp)
    Xp{i}(sum(X{i})==0)=NaN;
    Xp{i}(sum(X{i})==Inf)=NaN;
    Xp{i}(isnan(sum(X{i})))=NaN;
    temp=find(Xp{i}==1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extra dimentions
    X{i}=[X{i};zeros(length(ii),size(X{i},2))];
    for k=1:length(ii)
        X{i}(num_cost_terms+k,:)=X{i}(ii(k),:).*X{i}(jj(k),:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
[~,m]=size(X);

display(['Perceptron error-rate: ', num2str(nnz(w_perceptron*X<0)/m)])
display(['SVM error-rate: ', num2str(nnz(w_svm*X<0)/m)])
