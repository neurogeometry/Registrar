% This function calculates 5 cost components for merging branch pairs. These
% components are: tip-to-tip distance Dist, cosine of tip orientations, branch offset,
% average intensity, and variation of intensity along the length of the optimally connecting path.
% MaxDistance is the maximum distance between branches for which merger will be considered
% trim is the size of the region near the xy faces of the stack where branch tips are not merged.
% This function also prompts the user to perform semi-automated branch merger.
% The result is used for perceptron learning to determine the coefficients
% in the cost function.

function InitialMergingCostsStr = InitialMergingCosts(Orig,AMlbl,r,R,Parameters)

sizeIm=size(Orig);
L=unique(AMlbl(AMlbl>0));

% step back stepsback steps to determine tip1V tip2V vertexes of branch endings
[tip1V tip2V]=StepBack(AMlbl,Parameters.BranchMerger.StepsBack*Parameters.Optimization.PointsPerVoxel);

tip1=zeros(1,length(L));
tip2=zeros(1,length(L));
for i=1:length(L)
    if ~isempty(tip1V{L(i)})
        tip1(L(i))=tip1V{L(i)}(1);
        tip2(L(i))=tip2V{L(i)}(1);
    end
end

% branch end orientations
n1=zeros(length(L),3);
n2=zeros(length(L),3);
for i=1:length(L)
    if ~isempty(tip1V{L(i)})
        n1(i,:)=r(tip1V{L(i)}(end),:)-r(tip1V{L(i)}(1),:);
        n1(i,:)=n1(i,:)./sum(n1(i,:).^2)^0.5;
        
        n2(i,:)=r(tip2V{L(i)}(end),:)-r(tip2V{L(i)}(1),:);
        n2(i,:)=n2(i,:)./sum(n2(i,:).^2)^0.5;
    end
end

% 4 distances and cosines are calculated: tip1-tip1, tip1-tip2, tip2-tip1, tip2-tip2
Dist=inf(length(L),length(L),4); % Geometric distance. Dist=inf is the most unfavorable situation for merging
Overrun=inf(length(L),length(L),4); % Signed distance between branch tips in the branch plane. Positive for crossing branches, zero otherwise. Overrun=inf is the most unfavorable situation for merging
Offset=inf(length(L),length(L),4); % Out of branch plane component of distance between branch tips. Offset=inf is the most unfavorable situation for merging
Cos=ones(length(L),length(L),4); % Cos=1 is the most unfavorable situation for merging
%MeanI=zeros(length(L),length(L),4); % Average intensity along the snake. MeanI=0 is the most unfavorable situation for merging
%CVI=inf(length(L),length(L),4); % Coefficient of variation of intensity along the snake. CVI=inf is the most unfavorable situation for merging

for i=1:length(L)
    if ~isempty(tip1V{L(i)})
        for j=(i+1):length(L)
            if ~isempty(tip1V{L(j)})
                n_outofplane=cross_AS(n1(i,:),n1(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
                %n_perpendicular=n1(i,:)+n1(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
                n_parallel=n1(i,:)-n1(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
                
                r12=r(tip1V{L(i)}(1),:)-r(tip1V{L(j)}(1),:);
                %r_outofplane=abs(sum(r12.*n_outofplane)); % out of plane component of r12
                %r_perpendicular=abs(sum(r12.*n_perpendicular)); % perpendicular component of r12
                %r_parallel=sum(r12.*n_parallel); % parallel component of r12
                Dist(i,j,1)=sum(r12.^2).^0.5;
                Cos(i,j,1)=sum(n1(i,:).*n1(j,:));
                Overrun(i,j,1)=-sum(r12.*n_parallel);
                Offset(i,j,1)=abs(sum(r12.*n_outofplane));
                
                n_outofplane=cross_AS(n1(i,:),n2(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
                %n_perpendicular=n1(i,:)+n2(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
                n_parallel=n1(i,:)-n2(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
                
                r12=r(tip1V{L(i)}(1),:)-r(tip2V{L(j)}(1),:);
                %r_outofplane=abs(sum(r12.*n_outofplane));
                %r_perpendicular=abs(sum(r12.*n_perpendicular));
                %r_parallel=sum(r12.*n_parallel);
                Dist(i,j,2)=sum(r12.^2).^0.5;
                Cos(i,j,2)=sum(n1(i,:).*n2(j,:));
                Overrun(i,j,2)=-sum(r12.*n_parallel);
                Offset(i,j,2)=abs(sum(r12.*n_outofplane));
                
                n_outofplane=cross_AS(n2(i,:),n1(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
                %n_perpendicular=n2(i,:)+n1(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
                n_parallel=n2(i,:)-n1(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
                
                r12=r(tip2V{L(i)}(1),:)-r(tip1V{L(j)}(1),:);
                %r_outofplane=abs(sum(r12.*n_outofplane));
                %r_perpendicular=abs(sum(r12.*n_perpendicular));
                %r_parallel=sum(r12.*n_parallel);
                Dist(i,j,3)=sum(r12.^2).^0.5;
                Cos(i,j,3)=sum(n2(i,:).*n1(j,:));
                Overrun(i,j,3)=-sum(r12.*n_parallel);
                Offset(i,j,3)=abs(sum(r12.*n_outofplane));
                
                n_outofplane=cross_AS(n2(i,:),n2(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
                %n_perpendicular=n2(i,:)+n2(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
                n_parallel=n2(i,:)-n2(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
                
                r12=r(tip2V{L(i)}(1),:)-r(tip2V{L(j)}(1),:);
                %r_outofplane=abs(sum(r12.*n_outofplane));
                %r_perpendicular=abs(sum(r12.*n_perpendicular));
                %r_parallel=sum(r12.*n_parallel);
                Dist(i,j,4)=sum(r12.^2).^0.5;
                Cos(i,j,4)=sum(n2(i,:).*n2(j,:));
                Overrun(i,j,4)=-sum(r12.*n_parallel);
                Offset(i,j,4)=abs(sum(r12.*n_outofplane));
            end
        end
    end
end

% Overrun is positive for crossing branches and zero otherwise.
Overrun(Overrun<0)=0;

% distant branch tips or branch tips near the faces of the stack will not be mearged
ind=tip1(L); ind(ind==0)=[];
edge_tip1=(r(ind,1)<=Parameters.BranchMerger.TrimTrace | r(ind,1)>sizeIm(1)-Parameters.BranchMerger.TrimTrace | r(ind,2)<=Parameters.BranchMerger.TrimTrace | r(ind,2)>sizeIm(2)-Parameters.BranchMerger.TrimTrace);
ind=tip2(L); ind(ind==0)=[];
edge_tip2=(r(ind,1)<=Parameters.BranchMerger.TrimTrace | r(ind,1)>sizeIm(1)-Parameters.BranchMerger.TrimTrace | r(ind,2)<=Parameters.BranchMerger.TrimTrace | r(ind,2)>sizeIm(2)-Parameters.BranchMerger.TrimTrace);
Dist(edge_tip1,:,1)=inf; Dist(:,edge_tip1,1)=inf;
Dist(edge_tip1,:,2)=inf; Dist(:,edge_tip2,2)=inf;
Dist(edge_tip2,:,3)=inf; Dist(:,edge_tip1,3)=inf;
Dist(edge_tip2,:,4)=inf; Dist(:,edge_tip2,4)=inf;

InitialMergingCostsStr.Dist=Dist;
InitialMergingCostsStr.Overrun=Overrun;
InitialMergingCostsStr.Offset=Offset;
InitialMergingCostsStr.Cos=Cos;
InitialMergingCostsStr.n1=n1;
InitialMergingCostsStr.n2=n2;

% ind_pass=find(Dist(:)<Parameters.MaxDistance);
% [i_pass,j_pass,k_pass]=ind2sub(size(Dist),ind_pass);
%
% AM_merge=sparse([]);
% r_merge=zeros(numel(cell2mat(tip1V(i_pass)'))+numel(cell2mat(tip2V(j_pass)')),3);
% R_merge=zeros(numel(cell2mat(tip1V(i_pass)'))+numel(cell2mat(tip2V(j_pass)')),1);
% count=0;
% for i=1:length(k_pass)
%     if k_pass(i)==1
%         tempV1=tip1V{L(i_pass(i))};
%         tempV2=tip1V{L(j_pass(i))};
%     elseif k_pass(i)==2
%         tempV1=tip1V{L(i_pass(i))};
%         tempV2=tip2V{L(j_pass(i))};
%     elseif k_pass(i)==3
%         tempV1=tip2V{L(i_pass(i))};
%         tempV2=tip1V{L(j_pass(i))};
%     elseif k_pass(i)==4
%         tempV1=tip2V{L(i_pass(i))};
%         tempV2=tip2V{L(j_pass(i))};
%     end
%
%     N1=length(tempV1);
%     N2=length(tempV2);
%     r_merge(count+1:count+N1+N2,:)=r([tempV1,tempV2],:);
%     R_merge(count+1:count+N1+N2)=R([tempV1,tempV2]);
%     AM_merge(end+1:end+N1,end+1:end+N1)=triu(ones(N1,N1),1)-triu(ones(N1,N1),2);
%     AM_merge(end+1:end+N2,end+1:end+N2)=triu(ones(N2,N2),1)-triu(ones(N2,N2),2);
%     AM_merge(count+1,count+N1+1)=1;
%     count=count+N1+N2;
% end

% Intensity along the optimal connections
% [AMlbl_merge,r_merge,~,I_snake]=Snake_All_norm(Orig,AM_merge,r_merge,R_merge,0,0,0,Parameters.pointsperum,Parameters.Nstep,Parameters.alpha,Parameters.betta,Parameters.sig_ips,Parameters.sig_bps,Parameters.sig_tps,0);
% %[AMlbl_merge,~,~,I_snake]=Snake_rR(Orig,AM_merge,r_merge,R_merge,0,0,0,Parameters.pointsperum,Parameters.Nstep,Parameters.alpha,Parameters.betta,Parameters.sig_ips,Parameters.sig_bps,Parameters.sig_tps,1);
%
% L_new=unique(AMlbl_merge(AMlbl_merge>0));
% AMlbl_merge_direct=triu(AMlbl_merge);
% % Mean intensity and CV
% for i=1:length(L_new)
%     [e1 e2]=find(AMlbl_merge_direct==L_new(i));
%     Ibranch_segs=(I_snake(e1)+I_snake(e2))./2;
%     MeanI(ind_pass(L_new(i)))=mean(Ibranch_segs);
%     CVI(ind_pass(L_new(i)))=(mean(Ibranch_segs.^2)-mean(Ibranch_segs)^2)^0.5/mean(Ibranch_segs);
% end

%InitialMergingCostsStr.MeanI=MeanI;
%InitialMergingCostsStr.CVI=CVI;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % r_pre_end=zeros(length(k_pass),3);
% % r_post_end=zeros(length(k_pass),3);
% % r_pre_start=zeros(length(k_pass),3);
% % r_post_start=zeros(length(k_pass),3);
% % for i=1:length(k_pass)
% %     if k_pass(i)==1
% %         tempV1=tip1V{L(i_pass(i))};
% %         tempV2=tip1V{L(j_pass(i))};
% %     elseif k_pass(i)==2
% %         tempV1=tip1V{L(i_pass(i))};
% %         tempV2=tip2V{L(j_pass(i))};
% %     elseif k_pass(i)==3
% %         tempV1=tip2V{L(i_pass(i))};
% %         tempV2=tip1V{L(j_pass(i))};
% %     elseif k_pass(i)==4
% %         tempV1=tip2V{L(i_pass(i))};
% %         tempV2=tip2V{L(j_pass(i))};
% %     end
% %     r_pre_start(i,:)=r(tempV1(1),:);
% %     r_post_start(i,:)=r(tempV2(1),:);
% %     r_pre_end(i,:)=r(tempV1(end),:);
% %     r_post_end(i,:)=r(tempV2(end),:);
% % end
% %
% % % Plot the meargers
% % X=[Dist(ind_pass(L_new)),Overrun(ind_pass(L_new)),Offset(ind_pass(L_new)),Cos(ind_pass(L_new)),MeanI(ind_pass(L_new)),CVI(ind_pass(L_new))]';
% % disp(['Numbe of putative mergers is ',num2str(length(ind_pass))])
% % path='C:\Armen\DIADEM\Neuron Tracer V11\Parameter Files\Initial_XXp_data\';
% % figure(20)
% % imshow(max(Orig,[],3));
% % hold all; drawnow
% % PlotAM(AMlbl,r);
% % for i=1:length(L_new),i
% %     h_pre_start=plot3(r_pre_start(L_new(i),2),r_pre_start(L_new(i),1),r_pre_start(L_new(i),3),'go');
% %     h_pre_end=plot3(r_pre_end(L_new(i),2),r_pre_end(L_new(i),1),r_pre_end(L_new(i),3),'r*');
% %     h_post_start=plot3(r_post_start(L_new(i),2),r_post_start(L_new(i),1),r_post_start(L_new(i),3),'go');
% %     h_post_end=plot3(r_post_end(L_new(i),2),r_post_end(L_new(i),1),r_post_end(L_new(i),3),'r*');
% %     disp(['D=',num2str(Dist(ind_pass(L_new(i)))), '; Overrun=',num2str(Overrun(ind_pass(L_new(i)))),'; Offset=',num2str(Offset(ind_pass(L_new(i)))),'; Cos=',num2str(Cos(ind_pass(L_new(i))))])
% %     AM_merge_temp=sparse(size(AMlbl_merge,1),size(AMlbl_merge,2));
% %     AM_merge_temp(AMlbl_merge==L_new(i))=1;
% %     h3=PlotAM_h(AM_merge_temp,r_merge);
% %     user_input = input('Merge Type (2,3,4): ');
% %     MergerType(i)=user_input;
% %     user_input = input('Merge? (y,n,?): ','s');
% %     if user_input=='y'
% %         Xp(i)=1;
% %     elseif user_input=='n'
% %         Xp(i)=-1;
% %     else
% %         Xp(i)=NaN;
% %     end
% %     %save([path,'L6_06'],'X','Xp','MergerType','Parameters')
% %     delete(h_pre_start,h_pre_end,h_post_start,h_post_end,cell2mat(h3))
% % end
% %
% % %w2=perceptron(X,Xp,MergerType,2,0.1)




