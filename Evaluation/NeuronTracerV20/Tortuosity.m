% This function calculates maximum tortuosity, L/R-1, for every tree contained in AMlbl.
% Tmax is calculated in the [Lmin Lmax] range of branch lengths. 

function Tmax=Tortuosity(AMlbl,r,Lmin,Lmax)

Labels=unique(AMlbl(AMlbl>0));
Tmax=zeros(1,max(Labels));

for i=1:length(Labels)
    [e1,e2]=find(AMlbl==Labels(i));
    e12=unique([e1;e2]);
    AM_tree=AMlbl(e12,e12);
    AM_tree(AM_tree~=Labels(i))=0;
    AM_tree=((AM_tree+AM_tree')>0);
    r_tree=r(e12,:);
    
    % disconnect tree branches
    bp=find(sum(AM_tree)>2);
    for j=1:length(bp)
        bp_neigh=find(AM_tree(:,bp(j)));
        r_tree=[r_tree;ones(length(bp_neigh)-1,1)*r_tree(bp(j),:)];
        a=length(AM_tree);
        AM_tree(a+length(bp_neigh)-1,a+length(bp_neigh)-1)=0;
        AM_tree(bp_neigh(2:end),bp(j))=0;
        AM_tree(bp(j),bp_neigh(2:end))=0;
        temp_ind=sub2ind(size(AM_tree),bp_neigh(2:end),a-1+[2:length(bp_neigh)]');
        AM_tree(temp_ind)=1;
        temp_ind=sub2ind(size(AM_tree),a-1+[2:length(bp_neigh)]',bp_neigh(2:end));
        AM_tree(temp_ind)=1;
    end
    
    cv=find(sum(AM_tree)==1,1,'first');
    
    CumL=zeros(1,size(AM_tree,1));
    r_temp=zeros(size(r_tree));
    r_temp(1,:)=r_tree(cv,:);
    for j=2:size(AM_tree,1)
        nv=find(AM_tree(cv,:));
        if ~isempty(nv)
            CumL(j)=CumL(j-1)+sum((r_tree(nv,:)-r_tree(cv,:)).^2,2).^0.5;
        else
            CumL(j)=0;
            nv=find(sum(AM_tree)==1,1,'first');
        end
        r_temp(j,:)=r_tree(nv,:);
        AM_tree(cv,nv)=0;
        AM_tree(nv,cv)=0;
        cv=nv;
    end
    
    br_ind=[find(CumL==0),length(CumL)+1];
    for k=1:length(br_ind)-1
        temp=triu(ones(br_ind(k+1)-br_ind(k),br_ind(k+1)-br_ind(k)),1);
        ind=find(temp);
        [ii,jj]=ind2sub(size(temp),ind);
        L=zeros(1,length(ind));
        R=zeros(1,length(ind));
        for j=1:length(ind)
            L(j)=abs(CumL(br_ind(k)-1+ii(j))-CumL(br_ind(k)-1+jj(j)));
            R(j)=sum((r_temp(br_ind(k)-1+ii(j),:)-r_temp(br_ind(k)-1+jj(j),:)).^2,2).^0.5;
        end
        ind2=(L>=Lmin & L<=Lmax);
        if nnz(ind2)>0
            Tmax(Labels(i))=max(Tmax(Labels(i)),max(L(L>=Lmin & L<=Lmax)./R(L>=Lmin & L<=Lmax)-1));
        end
    end
end





