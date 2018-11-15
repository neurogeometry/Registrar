% This function generates the Clusters structure necessary for interactive training. 
% The input is initial trace, i.e. disconnected unbranched neurites.
% The structure contains information related to the 20 most likely merger scenarios
% within each cluster. It includes end point positions, labels, and ids for each cluster,
% as well as merger AM matrices and cost components for each scenario.
% Information about the lowest cost merger (based on the initial cost) within each cluster 
% is contained in the field best_merger (1 vs -1).

function ClustersStr=GenerateClustersStr(Im,AMlbl,r,R,InitialMergingCostsStr,MeanI0,Parameters)

ClustersStr=struct('end_point_ind',{},'end_point_label',{},'end_point_r',{},'scenarios',{},'cost_components',{},'best_merger',{});

Dist=InitialMergingCostsStr.Dist;
Cos=InitialMergingCostsStr.Cos;
Offset=InitialMergingCostsStr.Offset;
MeanI=InitialMergingCostsStr.MeanI;
CVI=InitialMergingCostsStr.CVI;

L=unique(AMlbl(AMlbl>0));
AM=spones(AMlbl);

% find branch end points
tip1V0=zeros(1,length(L));
tip2V0=zeros(1,length(L));
Blengths=zeros(1,length(L));
for i=1:length(L)
    [v1 ~]=find(AMlbl==L(i));
    tip12=sort(v1(sum(AM(:,v1))>=3 | sum(AM(:,v1))==1));
    tip1V0(L(i))=tip12(1);
    tip2V0(L(i))=tip12(2);
    Blengths(L(i))=length(unique(v1));
end

% step back stepsback steps to determine tip1V tip2V vertexes of branch endings 
[tip1V tip2V]=StepBack(AMlbl,Parameters.stepsback*Parameters.pointsperum);

End_dist=[[Dist(:,:,1),Dist(:,:,2)];[Dist(:,:,3),Dist(:,:,4)]];
End_cos=[[Cos(:,:,1),Cos(:,:,2)];[Cos(:,:,3),Cos(:,:,4)]];
End_Offset=[[Offset(:,:,1),Offset(:,:,2)];[Offset(:,:,3),Offset(:,:,4)]];
End_meanI=[[MeanI(:,:,1),MeanI(:,:,2)];[MeanI(:,:,3),MeanI(:,:,4)]];
End_CVI=[[CVI(:,:,1),CVI(:,:,2)];[CVI(:,:,3),CVI(:,:,4)]];

% Cluster the branch end-points by dist
End_dist_temp=zeros(size(End_dist)); End_dist_temp(End_dist<=Parameters.Dist_thr)=1;
[ClusterN,ClusterV] = FindClustersAM(spones(End_dist_temp+End_dist_temp'));

% Sort clusters by size: 2,4,3,5,6,7,...
ClusterSize=zeros(1,ClusterN);
for i=1:ClusterN
    ClusterSize(i)=length(ClusterV{i});
end
[ClusterSize,Cluster_ind]=sort(ClusterSize);
ClusterV=ClusterV(Cluster_ind);
Cluster_ind=[find(ClusterSize==2),find(ClusterSize==4),find(ClusterSize==3),find(ClusterSize>4)];
ClusterV=ClusterV(Cluster_ind);
ClusterN=length(ClusterV);

% Calculate the costs of all possible mergers
for i=1:ClusterN  
    disp([num2str(i/ClusterN),', ',num2str(length(ClusterV{i})),'-way merger'])
    TipLabels=ClusterV{i};
    TipLabels(TipLabels>size(Dist,1))=TipLabels(TipLabels>size(Dist,1))-size(Dist,1);
    ClustersStr(i).end_point_label=TipLabels;
    [Costs MergerAM]=NwayMerger(TipLabels,End_dist(ClusterV{i},ClusterV{i}),End_cos(ClusterV{i},ClusterV{i}),End_meanI(ClusterV{i},ClusterV{i}),End_CVI(ClusterV{i},ClusterV{i}),MeanI0,Parameters);
    
    % calculation of true tortuosity, Tmax, total intensity and length, Itotal and Ltotal, for each merger
    Tmax=inf(min(length(Costs),20),1); %inf(length(Costs),1); %!!!!!!!!!!!!!!!!!!!!!!!!!
    Itotal=zeros(size(Tmax));
    Ltotal=zeros(size(Tmax));
    Imean=zeros(size(Tmax));
    CV=zeros(size(Tmax));
    Kmax=zeros(size(Tmax));
    Kmean=zeros(size(Tmax));
    Free_tips=zeros(size(Tmax));
    ppL=ClusterV{i};
    
    ClustersStr(i).scenarios=MergerAM(:,:,1:length(Tmax));
    ClustersStr(i).best_merger=-ones(1,length(Tmax));
    %ClustersStr(i).best_merger(1)=1; % scenarios are ordered according to the initial cost. Branches of the initial trace must be disconnected 
    
    b_length=zeros(1,length(TipLabels));
    r_connect=zeros(numel(cell2mat(tip1V(TipLabels)')),3);
    count=0;
    AM_connect_temp=zeros(size(r_connect,1),size(r_connect,1));
    for k=1:length(TipLabels)
        if ppL(k)<=size(Dist,1)
            tempV=tip1V{TipLabels(k)};
            ClustersStr(i).end_point_ind(k)=tip1V0(TipLabels(k));
        else
            tempV=tip2V{TipLabels(k)};
            ClustersStr(i).end_point_ind(k)=tip2V0(TipLabels(k));
        end
        b_length(k)=length(tempV);
        r_connect(count+1:count+b_length(k),:)=r(tempV,:);

        old_label_ind=find(TipLabels(1:k-1)==TipLabels(k),1);
        N_same_Vs=b_length(old_label_ind)+b_length(k)-Blengths(TipLabels(k));
        if isempty(old_label_ind) || (~isempty(old_label_ind) && N_same_Vs<=0)
            AM_connect_temp(count+1:count+b_length(k),count+1:count+b_length(k))=triu(ones(b_length(k),b_length(k)),1)-triu(ones(b_length(k),b_length(k)),2);
        else % for short branches tip1V and tip2V are overlapping and have to be merged  
            AM_connect_temp(count+1:count+b_length(k)-N_same_Vs+1,count+1:count+b_length(k)-N_same_Vs+1)=triu(ones(b_length(k)-N_same_Vs+1,b_length(k)-N_same_Vs+1),1)-triu(ones(b_length(k)-N_same_Vs+1,b_length(k)-N_same_Vs+1),2);
            AM_connect_temp(sum(b_length(1:old_label_ind)),count+b_length(k)-N_same_Vs+1)=1;
        end
        count=count+b_length(k);
    end
    ClustersStr(i).end_point_r=r(ClustersStr(i).end_point_ind,:);
    cumb_length=1+cumsum([0,b_length(1:end-1)]);

      % This is 30% faster. Not finished.
%     [iii,jjj,kkk]=find(AM_connect_temp);
%     iii=iii*ones(1,length(Tmax))+ones(length(iii),1)*(length(AM_connect_temp).*[0:length(Tmax)-1]);
%     jjj=jjj*ones(1,length(Tmax))+ones(length(jjj),1)*(length(AM_connect_temp).*[0:length(Tmax)-1]);
%     kkk=kkk*ones(1,length(Tmax))+ones(length(kkk),1)*[0:length(Tmax)-1];
%     AM_connect_all=sparse(iii(:),jjj(:),kkk(:),length(AM_connect_temp)*length(Tmax),length(AM_connect_temp)*length(Tmax));
%     [AM_connect_all,r_connect_all,I_snake_all]=Snake_All_norm(Im,AM_connect_all,repmat(r_connect,length(Tmax),1),1,0,0,Parameters.pointsperum,Parameters.Nstep,Parameters.alpha,Parameters.betta,Parameters.sig_ips,Parameters.sig_bps,Parameters.sig_tps,0);
         
    for j=1:length(Tmax)      
        AM_connect=AM_connect_temp;
        [ii,jj]=find(MergerAM(:,:,j));
        AM_connect(sub2ind(size(AM_connect),cumb_length(ii),cumb_length(jj)))=1;
        
        if Loops(AM_connect)==0
            [AM_connect1,r_connect1,I_snake]=Snake_All_norm(Im,AM_connect,r_connect,1,0,0,Parameters.pointsperum,Parameters.Nstep,Parameters.alpha,Parameters.betta,Parameters.sig_ips,Parameters.sig_bps,Parameters.sig_tps,0);
            AM_connect1=(AM_connect1>0);
            AM_connect1 = LabelTreesAM(AM_connect1);
            
            Tmax(j)=max(Tortuosity(AM_connect1,r_connect1,0,50));
            
            [it,jt]=find(AM_connect1);
            lll=sum((r_connect1(it,:)-r_connect1(jt,:)).^2,2).^0.5;
            Ltotal(j)=sum(lll)/2;
            Itotal(j)=sum(lll'.*(I_snake(it)+I_snake(jt)))/4;
            Imean(j)=mean(I_snake);
            CV(j)=std(I_snake)/mean(I_snake);
            Free_tips(j)=nnz(sum(MergerAM(:,:,j)+MergerAM(:,:,j)')==0);
            
            [~,Ktotal_temp,Kmean_temp]=Curvature(AM_connect1,r_connect1,1);
            Kmax(j)=max(Kmean_temp); % maximum curvature of connecting trees
            Kmean(j)=sum(Ktotal_temp)/Ltotal(j); % average curvature of all connecting trees
        end   
    end
    ClustersStr(i).cost_components=[Ltotal,Itotal,Imean,Tmax,CV,Free_tips,Kmax,Kmean]';
    
    % remove non-numeric costs scenarios
    rem=~isfinite(sum(ClustersStr(i).cost_components,1));
    ClustersStr(i).scenarios(:,:,rem)=[];
    ClustersStr(i).cost_components(:,rem)=[];
    ClustersStr(i).best_merger(rem)=[];
end

% remove clusters containing no scenarios (> 12 end-points)
rem=zeros(1,ClusterN);
for i=1:ClusterN
    rem(i)=isempty(ClustersStr(i).scenarios);
end
ClustersStr(rem==1)=[];
