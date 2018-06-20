% This function calculates 4 cost components for merging branch pairs. These 
% components are: tip-to-tip distance Dist, cosine of tip orientations,
% average intensity and variation of intensity along the length of the optimally connecting segment.
% D_thr and Cos_thr are the maximum distance between branches and maximum cosine
% for which merger will be considered
% This function also prompts the user to perform semi-automated branch merger.
% The result is used for perceptron learning to determine the coefficients
% in the cost function.
function [Dist,Cos,MeanI,FanoI] = MergingCosts_Sparse(Im,sizeIm,AMlbl,r,D_thr,stepsback,pointsperum,Nstep,alpha,betta)

% pointsperum=1;%/2;
% stepsback=50; % number of segments steped back to determine branch tip orientations and strech snakes betwen branches
AM=spones(AMlbl);
L=unique(AMlbl(AMlbl>0));

% tip1 tip2 are sorted (by vertex index) combination of branch tips
tip1=zeros(length(L),1);
tip2=zeros(length(L),1);
for i=1:length(L)   
    [v1 v2]=find(AMlbl==L(i));
    temp=sort(v1(sum(AM(:,v1))>=3 | sum(AM(:,v1))==1));
    tip1(L(i))=temp(1); % end 1
    tip2(L(i))=temp(2);
end
% rem_ind=(tip1(L)==tip2(L));
% tip1(L(rem_ind))=[];
% tip2(L(rem_ind))=[];
% L(rem_ind)=[];

% step back stepsback steps to determine orientations of branch endings at tip1 tip2
n1=zeros(length(L),3);
n2=zeros(length(L),3);
r1_back=zeros(length(L),3);
r2_back=zeros(length(L),3);
for i=1:length(L) 
    tempAM=(AMlbl==L(i));
    st=stepsback;
    if nnz(tempAM)/2<stepsback
        st=nnz(tempAM)/2;
    end
    
    v=tip1(L(i));
    if tip1(L(i))~=tip2(L(i))
        for f=1:st
            v_new=find(tempAM(v,:));
            tempAM(v,v_new)=0;
            tempAM(v_new,v)=0;
            v=v_new;
        end
    else
        v=find(tempAM(v,:));
        v=v(1);
    end
    %line([r(tip1(L(i)),2) r(v,2)],[r(tip1(L(i)),1) r(v,1)],[r(tip1(L(i)),3) r(v,3)],'Color','m')
    n1(i,:)=r(v,:)-r(tip1(L(i)),:);
    n1(i,:)=n1(i,:)./sum(n1(i,:).^2)^0.5;
    r1_back(i,:)=r(v,:);
    
    tempAM=(AMlbl==L(i));
    v=tip2(L(i));
    if tip1(L(i))~=tip2(L(i))
        for f=1:st
            v_new=find(tempAM(v,:));
            tempAM(v,v_new)=0;
            tempAM(v_new,v)=0;
            v=v_new;
        end
    else
        v=find(tempAM(v,:));
        v=v(2);
    end
    %line([r(tip2(L(i)),2) r(v,2)],[r(tip2(L(i)),1) r(v,1)],[r(tip2(L(i)),3) r(v,3)],'Color','m')
    n2(i,:)=r(v,:)-r(tip2(L(i)),:);
    n2(i,:)=n2(i,:)./sum(n2(i,:).^2)^0.5;
    r2_back(i,:)=r(v,:);
end

% 4 distances and cosines are calculated: tip1-tip1, tip1-tip2, tip2-tip1, tip2-tip2
Dist=inf(length(L),length(L),4); % Dist=inf is the most unfavorable situation for merging
Cos=ones(length(L),length(L),4); % Cos=1 is the most unfavorable situation for merging
MeanI=zeros(length(L),length(L),4); % MeanI=0 is the most unfavorable situation for merging
FanoI=inf(length(L),length(L),4); % FanoI=inf is the most unfavorable situation for merging

for i=1:length(L)
    for j=(i+1):length(L)
        % Distances
        Dist(i,j,1)=sum((r(tip1(L(i)),:)-r(tip1(L(j)),:)).^2).^0.5;
        Dist(i,j,2)=sum((r(tip1(L(i)),:)-r(tip2(L(j)),:)).^2).^0.5;
        Dist(i,j,3)=sum((r(tip2(L(i)),:)-r(tip1(L(j)),:)).^2).^0.5;
        Dist(i,j,4)=sum((r(tip2(L(i)),:)-r(tip2(L(j)),:)).^2).^0.5;
        
        % Cosines
        Cos(i,j,1)=sum(n1(i,:).*n1(j,:));
        Cos(i,j,2)=sum(n1(i,:).*n2(j,:));
        Cos(i,j,3)=sum(n2(i,:).*n1(j,:));
        Cos(i,j,4)=sum(n2(i,:).*n2(j,:));
    end
end

ind_pass=find(Dist(:)<D_thr);
[i_pass,j_pass,k_pass]=ind2sub(size(Dist),ind_pass);

AM_merge=sparse([]);
r_merge=[];
r_pre=zeros(length(k_pass),3);
r_post=zeros(length(k_pass),3);
for f=1:length(k_pass)
    if k_pass(f)==1
        r_pre(f,:)=r1_back(L(i_pass(f)),:); %r(tip1(L(i_pass(f))),:);
        r_post(f,:)=r1_back(L(j_pass(f)),:); %r(tip1(L(j_pass(f))),:);
    elseif k_pass(f)==2
        r_pre(f,:)=r1_back(L(i_pass(f)),:); %r(tip1(L(i_pass(f))),:);
        r_post(f,:)=r2_back(L(j_pass(f)),:); %r(tip2(L(j_pass(f))),:);
    elseif k_pass(f)==3
        r_pre(f,:)=r2_back(L(i_pass(f)),:); %r(tip2(L(i_pass(f))),:);
        r_post(f,:)=r1_back(L(j_pass(f)),:); %r(tip1(L(j_pass(f))),:);
    elseif k_pass(f)==4
        r_pre(f,:)=r2_back(L(i_pass(f)),:); %r(tip2(L(i_pass(f))),:);
        r_post(f,:)=r2_back(L(j_pass(f)),:); %r(tip2(L(j_pass(f))),:);
    end
    d=sum((r_pre(f,:)-r_post(f,:)).^2)^0.5;
    nsegs=ceil(d*pointsperum);
    if nsegs==0
        nsegs=1;
    end
    dr=[0:nsegs]'*(r_post(f,:)-r_pre(f,:))./nsegs;
    r_merge=[r_merge; [r_pre(f,1)+dr(:,1),r_pre(f,2)+dr(:,2),r_pre(f,3)+dr(:,3)]];
    AM_merge(end+1:end+nsegs+1,end+1:end+nsegs+1)=triu(ones(nsegs+1,nsegs+1),1)-triu(ones(nsegs+1,nsegs+1),2);
end

% Intensity along the optimal connections
[AMlbl_merge,temp,I_snake]=Snake_subd_AM2_Sparse(Im,sizeIm,AM_merge,r_merge,pointsperum,Nstep,alpha,betta);

AMlbl_merge_direct=triu(AMlbl_merge);
L_new=unique(AMlbl_merge_direct(AMlbl_merge_direct>0));
for i=1:length(L_new)
    [e1 e2]=find(AMlbl_merge_direct==L_new(i));
    Ibranch_segs=(I_snake(e1)+I_snake(e2))./2;
    MeanI(ind_pass(L_new(i)))=mean(Ibranch_segs);
    FanoI(ind_pass(L_new(i)))=(mean(Ibranch_segs.^2)-mean(Ibranch_segs)^2)/mean(Ibranch_segs);
end

