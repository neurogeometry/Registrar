% This function performs automated branch merger in the order of decreasing confidence.
% Cost for each merger is calculated beased on cost_components contained in
% ClustersStr and weights contained in the TrainingHistory structure. AM_merged is the 
% adjacency matrix for the merged branched structure. 
% Loops are avoided at every step of the merger. 

function [AMlbl_merged,r,R]=AutomatedMerger(AMlbl,r,R,ClustersStr,TrainingHistory,Parameters)

disp('Automated branch merger started.')

T=0.5; % temperature used to define the confidence measure
if ~isfield(TrainingHistory,'w')
    error('Incorrect Training History file format.')
elseif isempty(TrainingHistory(1).w)
    error('There is no training history.')
else
    Weights=TrainingHistory(1).w;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort ClustersStr scenarios according to confidence
Q=zeros(size(ClustersStr,2),3); %[confidence, cluster number, scenario number]
X_0={ClustersStr.cost_components};
count=0;
for i=1:size(ClustersStr,2) 
    Cost=Weights*X_0{i};
    confidence=exp((min(Cost)-Cost)./T)/sum(exp((min(Cost)-Cost)./T));
    temp=size(ClustersStr(i).scenarios,3);
    Q(count+1:count+temp,1)=confidence';
    Q(count+1:count+temp,2)=ones(temp,1).*i;
    Q(count+1:count+temp,3)=(1:temp)';
    count=count+temp;
end
Q=sortrows(Q,-1);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_Labels=max([ClustersStr.end_point_label]);
if isempty(N_Labels)
    N_Labels=0;
end
AMmergerL=sparse(N_Labels,N_Labels);
while ~isempty(Q)
    % merge scenarios based on confidence
    [pre,post]=find(ClustersStr(Q(1,2)).scenarios(:,:,Q(1,3)));
    if isempty(pre)
        Q(Q(:,2)==Q(1,2),:)=[];
    else
        preL=ClustersStr(Q(1,2)).end_point_label(pre);
        postL=ClustersStr(Q(1,2)).end_point_label(post);
        
        VpreL=ClustersStr(Q(1,2)).end_point_ind(pre);
        VpostL=ClustersStr(Q(1,2)).end_point_ind(post);
        
        temp_ind1=sub2ind(size(AMmergerL),preL,postL);
        temp_ind2=sub2ind(size(AMmergerL),postL,preL);
        if length(unique(temp_ind1))<length(temp_ind1) || nnz(AMmergerL(temp_ind1))>0 || nnz(AMmergerL(temp_ind2))>0
            Q(1,:)=[];
        else
            AMtemp=AMmergerL;
            AMtemp(temp_ind1)=1;
            if nnz(AMtemp+AMtemp')==2*nnz(AMtemp) && (Parameters.BranchMerger.AllowLoops==1 || Loops(AMtemp)==0)
                AMlbl_temp=AMlbl;
                temp_ind=(VpreL~=VpostL);
                AMlbl_temp(sub2ind(size(AMlbl),VpreL(temp_ind),VpostL(temp_ind)))=preL(temp_ind);
                
                [preL_temp,postL_temp]=find(AMtemp);
                for ii=1:length(preL_temp)
                    AMlbl_temp(AMlbl_temp==postL_temp(ii))=preL_temp(ii);
                    preL_temp(preL_temp==postL_temp(ii))=preL_temp(ii);
                    postL_temp(postL_temp==postL_temp(ii))=preL_temp(ii);
                end
                
                AMmergerL=AMtemp;
                AMlbl(sub2ind(size(AMlbl),VpreL(temp_ind),VpostL(temp_ind)))=preL(temp_ind);
                uVpreL=unique(VpreL(temp_ind));
                uVpostL=VpostL(temp_ind);
                for j=1:length(uVpreL)
                    tempind=uVpostL(VpreL(temp_ind)==uVpreL(j));
                    r(uVpreL(j),:)=mean([r(uVpreL(j),:);r(tempind,:)]);
                    R(uVpreL(j))=max([R(uVpreL(j));R(tempind)]);
                end
                Q(Q(:,2)==Q(1,2),:)=[];
            else
                Q(1,:)=[];
            end
        end
    end
end

% Create undirected AMlbl
AMlblp=AMlbl';
ind=(AMlbl-AMlblp)<0;
AMlbl(ind)=AMlblp(ind);

% Relabel AMlbl
[preL,postL]=find(AMmergerL);
AMlbl_merged=AMlbl;
for i=1:length(preL)
    AMlbl_merged(AMlbl_merged==postL(i))=preL(i);
    preL(preL==postL(i))=preL(i);
    postL(postL==postL(i))=preL(i);
end

% Reduce the lablels
L=sort(unique(AMlbl_merged(AMlbl_merged>0)));
for i=1:length(L)
    AMlbl_merged(AMlbl_merged==L(i))=i;
end

disp('Automated branch merger is complete.')

