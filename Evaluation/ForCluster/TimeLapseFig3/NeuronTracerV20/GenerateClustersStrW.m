% This function generates the Clusters structure necessary for interactive training. 
% The input is initial trace, i.e. disconnected unbranched neurites.
% The structure contains information related to the Parameters.NBestScenarios most likely merger scenarios
% within each cluster. It includes end point positions, labels, and ids for each cluster,
% as well as merger AM matrices, cost components for each scenario, weights reflecting user confidence.
% Information about the lowest cost merger (based on the initial cost) within each cluster 
% is contained in the field best_merger (1 vs -1).

% FINAL COST COMPONENTS (end-points must be fixed during optimization)
% 1. Dist (combined distance, for connecting end-points)
% 2. Overrun (combined overrun, for connecting end-points)
% 3. Offset (combined offset, for connecting end-points)
% 4. Cos2, 2-branch deviation from perfect branching angles (combined L1 norm)
% 5. Cos3, 3-branch deviation from perfect branching angles (combined L1 norm)
% 6. Cos4, 4-branch deviation from perfect branching angles (combined L1 norm)
% -. Plan, planarity of branching
% -. Zdiv, combined standard deviation for all trees
% 9. Icv, coefficient of variation in intensity (mean of CVs for every tree)
% 10. Rcv, coefficient of variation in radius (mean of CVs for every tree)
% 11. N1, number of free end-points
% -. N2, number of 2-point mergers
% -. N3, number of byfurcations
% -. N4, number of tryfurcations

function [AMlbl,r,R,ClustersStr]=GenerateClustersStrW(Orig,AMlbl,r,R,Parameters,varargin)

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

disp('Calculation of cost components started.')

% Disconnect all branches at branch points (if not already disconnected)
[AMlbl,r,R] = Disconnect_Branches(AMlbl,r,R);

% calculation of intensities
[~,~,~,I_F]=Optimize_Trace_R(Orig,AMlbl,r,R,0,1,0,Parameters.Optimization.PointsPerVoxel,1,Parameters.Optimization.TraceStiffness,Parameters.Optimization.RadiusStiffness,0,0,0,myClass);

if myClass.getFlag()==1
    return;
end

ClustersStr=struct('end_point_ind',{},'end_point_label',{},'end_point_r',{},'scenarios',{},'cost_components',{},'best_merger',{},'alpha',{},'w',{},'weights',{},'confidence',{});

L=unique(AMlbl(AMlbl>0));

Cos2best=-1; Cos3best=-0.5; Cos4best=0;

% step back stepsback steps to determine tip1V tip2V vertexes of branch endings 
[tip1V tip2V]=StepBack(AMlbl,Parameters.BranchMerger.StepsBack*Parameters.Optimization.PointsPerVoxel);

% find branch end points
tip1V0=zeros(1,length(L));
tip2V0=zeros(1,length(L));
Blengths=zeros(1,length(L));
for i=1:length(L)
    [v1 ~]=find(AMlbl==L(i));
    if ~isempty(tip1V{L(i)})
        tip1V0(L(i))=tip1V{L(i)}(1);
        tip2V0(L(i))=tip2V{L(i)}(1);
    end
    Blengths(L(i))=length(unique(v1));
end

InitialMergingCostsStr = InitialMergingCosts(Orig,AMlbl,r,R,Parameters);

End_dist=[[InitialMergingCostsStr.Dist(:,:,1),InitialMergingCostsStr.Dist(:,:,2)];[InitialMergingCostsStr.Dist(:,:,3),InitialMergingCostsStr.Dist(:,:,4)]];
End_overrun=[[InitialMergingCostsStr.Overrun(:,:,1),InitialMergingCostsStr.Overrun(:,:,2)];[InitialMergingCostsStr.Overrun(:,:,3),InitialMergingCostsStr.Overrun(:,:,4)]];
End_offset=[[InitialMergingCostsStr.Offset(:,:,1),InitialMergingCostsStr.Offset(:,:,2)];[InitialMergingCostsStr.Offset(:,:,3),InitialMergingCostsStr.Offset(:,:,4)]];
End_cos=[[InitialMergingCostsStr.Cos(:,:,1),InitialMergingCostsStr.Cos(:,:,2)];[InitialMergingCostsStr.Cos(:,:,3),InitialMergingCostsStr.Cos(:,:,4)]];

End_dist=min(End_dist,End_dist');
End_overrun=min(End_overrun,End_overrun');
End_offset=min(End_offset,End_offset');
End_cos=min(End_cos,End_cos');

% branch orientations
%End_direction=[InitialMergingCostsStr.n1;InitialMergingCostsStr.n2];

% Cluster the branch end-points by dist
End_dist_temp=sparse(size(End_dist,1),size(End_dist,2)); 
End_dist_temp(End_dist<=Parameters.BranchMerger.MinClusterDist)=1;
%%%%%%%%%%%%% (???) Can in addition threshold clusters based on End_overrun and End_offset %%%%%%%%%%%%% 
[ClusterN,ClusterV] = FindClustersAM(End_dist_temp);

%%%%%%%%%%%%% (2) Can subdivide big clusters based on End_dist, End_overrun, and End_offset %%%%%%%%%%%%%
%%%%%%%%%%%%% (3) Can subdivide big clusters based on Initial Cost %%%%%%%%%%%%%

% Sort clusters by size: 2,3,4,5,6,7,...
ClusterSize=zeros(1,ClusterN);
for i=1:ClusterN
    ClusterSize(i)=length(ClusterV{i});
end
[~,Cluster_ind]=sort(ClusterSize);
ClusterV=ClusterV(Cluster_ind);
ClusterN=length(ClusterV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % training module
% figure(100)
% imshow(max(Orig,[],3),[0 max(Orig(:))]) 
% hold on 
% PlotAM(AMlbl, r)
% 
% figure(101)
% imshow(max(Orig,[],3),[0 max(Orig(:))]) 
% hold on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the costs of all possible mergers
for i=1:ClusterN
    disp([num2str(i/ClusterN),', ',num2str(length(ClusterV{i})),'-point merger'])
    TipLabels=ClusterV{i};
    TipLabels(TipLabels>length(L))=TipLabels(TipLabels>length(L))-length(L);
    ClustersStr(i).end_point_label=TipLabels;
    [Costs MergerAM]=NwayMerger(TipLabels,End_dist(ClusterV{i},ClusterV{i}),End_overrun(ClusterV{i},ClusterV{i}),End_offset(ClusterV{i},ClusterV{i}),End_cos(ClusterV{i},ClusterV{i}),Parameters);
    
    Nmergers=min(length(Costs),Parameters.BranchMerger.NBestScenarios); %length(Costs) %%%%%%%%%%%%%%%%%%%%%%%%%

    ClustersStr(i).scenarios=MergerAM(:,:,1:Nmergers);
    ClustersStr(i).best_merger=-ones(1,Nmergers);
    %ClustersStr(i).best_merger(1)=1; % scenarios are ordered according to the initial cost. Branches of the initial trace must be disconnected 
    ClustersStr(i).alpha=zeros(1,Nmergers);
    %ClustersStr(i).w=[];
    ClustersStr(i).weights=0;
    
    b_length=zeros(1,length(TipLabels));
    r_connect_temp=zeros(numel(cell2mat(tip1V(TipLabels)')),3);
    R_connect_temp=zeros(numel(cell2mat(tip1V(TipLabels)')),1);
    I_F_connect_temp=zeros(numel(cell2mat(tip1V(TipLabels)')),1);
    count=0;
    AM_connect_temp=zeros(size(r_connect_temp,1),size(r_connect_temp,1));
    for k=1:length(TipLabels)
        if myClass.getFlag()==1
            return;
        end

        if ClusterV{i}(k)<=length(L)
            tempV=tip1V{TipLabels(k)};
            ClustersStr(i).end_point_ind(k)=tip1V0(TipLabels(k));
        else
            tempV=tip2V{TipLabels(k)};
            ClustersStr(i).end_point_ind(k)=tip2V0(TipLabels(k));
        end
        b_length(k)=length(tempV);
        r_connect_temp(count+1:count+b_length(k),:)=r(tempV,:);
        R_connect_temp(count+1:count+b_length(k))=R(tempV);
        I_F_connect_temp(count+1:count+b_length(k))=I_F(tempV);
        
        old_label_ind=find(TipLabels(1:k-1)==TipLabels(k),1);
        N_same_Vs=b_length(old_label_ind)+b_length(k)-Blengths(TipLabels(k));
        if isempty(old_label_ind) || (~isempty(old_label_ind) && N_same_Vs<=0)
            AM_connect_temp(count+1:count+b_length(k),count+1:count+b_length(k))=triu(ones(b_length(k),b_length(k)),1)-triu(ones(b_length(k),b_length(k)),2);
        elseif b_length(k)>N_same_Vs % for short branches tip1V and tip2V are overlapping and have to be merged             
            AM_connect_temp(count+1:count+b_length(k)-N_same_Vs,count+1:count+b_length(k)-N_same_Vs)=triu(ones(b_length(k)-N_same_Vs,b_length(k)-N_same_Vs),1)-triu(ones(b_length(k)-N_same_Vs,b_length(k)-N_same_Vs),2);
            AM_connect_temp(sum(b_length(1:old_label_ind)),count+b_length(k)-N_same_Vs)=1;
        elseif b_length(k)==N_same_Vs % creating zero-length segment
            AM_connect_temp(sum(b_length(1:old_label_ind)),count+b_length(k)-N_same_Vs+1)=1;
        end
        count=count+b_length(k);
    end
    ClustersStr(i).end_point_r=r(ClustersStr(i).end_point_ind,:);
    cumb_length=1+cumsum([0,b_length(1:end-1)]);

    Dist=nan(1,Nmergers);
    Overrun=nan(1,Nmergers);
    Offset=nan(1,Nmergers);
    Cos2=nan(1,Nmergers);
    Cos3=nan(1,Nmergers);
    Cos4=nan(1,Nmergers);
    %Plan=nan(1,Nmergers);
    %Zdiv=nan(1,Nmergers);
    Icv=nan(1,Nmergers);
    Rcv=nan(1,Nmergers);
    N1=nan(1,Nmergers);
    %N2=nan(1,Nmergers);
    %N3=nan(1,Nmergers);
    %N4=nan(1,Nmergers);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % training module
%     user_input='n';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j=1:Nmergers
        if myClass.getFlag()==1
            return;
        end
        
        AM_connect=AM_connect_temp;
        [ii,jj]=find(MergerAM(:,:,j));
        AM_connect(sub2ind(size(AM_connect),cumb_length(ii),cumb_length(jj)))=1;
        AM_connect=spones(AM_connect+AM_connect');
        
        if Loops(AM_connect)==0              
            % CALCULATION OF FINAL COST COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%
            ind=sub2ind(size(MergerAM(:,:,j)),ii,jj);
            MergerAMsym=MergerAM(:,:,j)+MergerAM(:,:,j)';
            temp=sum(MergerAMsym);
            
            Dist(j)=0;
            Overrun(j)=0;
            Offset(j)=0;
            Cos2(j)=0;
            Cos3(j)=0;
            Cos4(j)=0;
            %Plan(j)=0;
                
            if ~isempty(ind)  
                temp_dist=End_dist(ClusterV{i},ClusterV{i});
                temp_overrun=End_overrun(ClusterV{i},ClusterV{i});
                temp_offset=End_offset(ClusterV{i},ClusterV{i});
                temp_cos=End_cos(ClusterV{i},ClusterV{i});
                
                Dist(j)=sum(temp_dist(ind));
                Overrun(j)=sum(temp_overrun(ind));
                Offset(j)=sum(temp_offset(ind));
                
                e1=find(temp==1);
                for s=1:length(e1)
                    e2=find(MergerAMsym(e1(s),:));
                    if temp(e2)==1
                        Cos2(j)=Cos2(j)+abs(temp_cos(e1(s),e2)-Cos2best)/2;
                    end
                end
                                
                %temp_direction=End_direction(ClusterV{i},:);
                e1=find(temp==2);
                for s=1:length(e1)
                    e23=find(MergerAMsym(e1(s),:));
                    tempCos=triu(abs(temp_cos([e1(s),e23],[e1(s),e23])-Cos3best),1);
                    Cos3(j)=Cos3(j)+sum(tempCos(:))/3;
                    %Cnn=temp_direction([e1(s),e23],:)'*temp_direction([e1(s),e23],:)./3;
                    %Plan(j)=Plan(j)+min(eig(Cnn));
                end
                e1=find(temp==3);
                for s=1:length(e1)
                    e234=find(MergerAMsym(e1(s),:));
                    tempCos=triu(abs((temp_cos([e1(s),e234],[e1(s),e234])-Cos2best).*(temp_cos([e1(s),e234],[e1(s),e234])-Cos4best)),1);
                    Cos4(j)=Cos4(j)+sum(tempCos(:))/6;
                    %Cnn=temp_direction([e1(s),e234],:)'*temp_direction([e1(s),e234],:)./4;
                    %Plan(j)=Plan(j)+min(eig(Cnn));
                end
            end
            
            AM_connect=LabelTreesAM(AM_connect);
            L_connect=unique(AM_connect(AM_connect>0));
            %Ztree_var=zeros(1,length(L_connect));
            Itree_cv=zeros(1,length(L_connect));
            Rtree_cv=zeros(1,length(L_connect));
            for s=1:length(L_connect)
                [e1 ~]=find(AM_connect==L_connect(s));
                e1=unique(e1);
                %Ztree_var(s)=var(r_connect_temp(e1,3));
                Itree_cv(s)=std(I_F_connect_temp(e1))/mean(I_F_connect_temp(e1));
                Rtree_cv(s)=std(R_connect_temp(e1))/mean(R_connect_temp(e1));
            end 
            %Zdiv(j)=(mean(Ztree_var))^0.5;
            Icv(j)=mean(Itree_cv);
            Rcv(j)=mean(Rtree_cv);    

            N1(j)=sum(temp==0);
            %N3(j)=sum(temp==2);
            %N4(j)=sum(temp==3);
            %N2(j)=sum(temp)/2-N3(j)*2-N4(j)*3;
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % training module
%             if ~strcmp(user_input,'y') && ~strcmp(user_input,'stop')
%                 display([Dist(j),Overrun(j),Offset(j),Cos2(j),Cos3(j),Cos4(j),Icv(j),Rcv(j),N1(j)])
%                 
%                 figure(101)
%                 h1=PlotAM_h(AM_connect, r_connect_temp);
%                 h2=plot3(r_connect_temp(cumb_length,2),r_connect_temp(cumb_length,1),r_connect_temp(cumb_length,3),'r*');
%                 user_input = input('Merge y/n?: ','s');
%                 if user_input=='y'
%                     ClustersStr(i).best_merger(j)=1;
%                     ClustersStr(i).weights=1;
%                 elseif user_input=='n'
%                     ClustersStr(i).best_merger(j)=-1;
%                 else
%                     ClustersStr(i).best_merger(j)=NaN;
%                 end
%                 delete(cell2mat(h1'),h2) 
%             end
% %             TrainingHistory = WSVMClassifier(ClustersStr,[],Parameters);
% %             save('C:\Armen\DIADEM\Neuron Tracer V12\Parameter Files\TrainingHistory_L6','TrainingHistory','Parameters')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        end   
    end
    
    %ClustersStr(i).cost_components=[Dist;Overrun;Offset;Cos2;Cos3;Cos4;Plan;Zdiv;Icv;Rcv;N1;N2;N3;N4];
    ClustersStr(i).cost_components=[Dist;Overrun;Offset;Cos2;Cos3;Cos4;Icv;Rcv;N1];
    
    % remove non-numeric costs scenarios
    rem=(~isfinite(sum(ClustersStr(i).cost_components,1)) | isnan(sum(ClustersStr(i).cost_components,1)));
    ClustersStr(i).scenarios(:,:,rem)=[];
    ClustersStr(i).cost_components(:,rem)=[];
    ClustersStr(i).best_merger(rem)=[];
    ClustersStr(i).alpha(rem)=[];
    
    if ~isempty(ClustersStr(i).scenarios)
        ClustersStr(i).confidence=1/size(ClustersStr(i).scenarios,3);
    end
end

% remove clusters containing no scenarios (> 12 end-points)
rem=zeros(1,ClusterN);
for i=1:ClusterN
    rem(i)=isempty(ClustersStr(i).scenarios);
end
ClustersStr(rem==1)=[];

disp('Calculation of cost components is complete.')
