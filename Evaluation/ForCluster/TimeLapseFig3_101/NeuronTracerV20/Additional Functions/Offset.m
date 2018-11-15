% This function calculates the distance (D), cosine (C) and offset (O)
% for every branch pair contained in AMlbl

function [D C O]=Offset(AMlbl,r,stepsback)

Labels=unique(AMlbl(AMlbl>0));
D=zeros(max(Labels),max(Labels),4);
C=zeros(max(Labels),max(Labels),4);
O=zeros(max(Labels),max(Labels),4);

% step back stepsback steps to determine tip1V tip2V vertexes of branch endings 
[tip1V tip2V]=StepBack(AMlbl,stepsback);

% branch end orientations
n1=zeros(length(Labels),3);
n2=zeros(length(Labels),3);
for i=1:length(Labels) 
    n1(i,:)=r(tip1V{Labels(i)}(end),:)-r(tip1V{Labels(i)}(1),:);
    n1(i,:)=n1(i,:)./sum(n1(i,:).^2)^0.5;
    
    n2(i,:)=r(tip2V{Labels(i)}(end),:)-r(tip2V{Labels(i)}(1),:);
    n2(i,:)=n2(i,:)./sum(n2(i,:).^2)^0.5;
end

for i=1:length(Labels)
    for j=(i+1):length(Labels)
        %n_outofplane=cross_AS(n1(i,:),n1(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
        %n_perpendicular=n1(i,:)+n1(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
        n_parallel=n1(i,:)-n1(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
        
        r12=r(tip1V{Labels(i)}(1),:)-r(tip1V{Labels(j)}(1),:);
        %r_outofplane=abs(sum(r12.*n_outofplane)); % out of plane component of r12
        %r_perpendicular=abs(sum(r12.*n_perpendicular)); % perpendicular component of r12
        %r_parallel=sum(r12.*n_parallel); % parallel component of r12
        D(i,j,1)=sum(r12.^2).^0.5;
        C(i,j,1)=sum(n1(i,:).*n1(j,:));
        O(i,j,1)=sum(r12.*n_parallel);
        
        %n_outofplane=cross_AS(n1(i,:),n2(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
        %n_perpendicular=n1(i,:)+n2(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
        n_parallel=n1(i,:)-n2(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
        
        r12=r(tip1V{Labels(i)}(1),:)-r(tip2V{Labels(j)}(1),:);
        %r_outofplane=abs(sum(r12.*n_outofplane));
        %r_perpendicular=abs(sum(r12.*n_perpendicular));
        %r_parallel=sum(r12.*n_parallel);
        D(i,j,2)=sum(r12.^2).^0.5;
        C(i,j,2)=sum(n1(i,:).*n2(j,:));
        O(i,j,2)=sum(r12.*n_parallel);
        
        %n_outofplane=cross_AS(n2(i,:),n1(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
        %n_perpendicular=n2(i,:)+n1(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
        n_parallel=n2(i,:)-n1(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
        
        r12=r(tip2V{Labels(i)}(1),:)-r(tip1V{Labels(j)}(1),:);
        %r_outofplane=abs(sum(r12.*n_outofplane));
        %r_perpendicular=abs(sum(r12.*n_perpendicular));
        %r_parallel=sum(r12.*n_parallel);
        D(i,j,3)=sum(r12.^2).^0.5;
        C(i,j,3)=sum(n2(i,:).*n1(j,:));
        O(i,j,3)=sum(r12.*n_parallel);
        
        %n_outofplane=cross_AS(n2(i,:),n2(j,:)); n_outofplane=n_outofplane./sum(n_outofplane.^2)^0.5;
        %n_perpendicular=n2(i,:)+n2(j,:); n_perpendicular=n_perpendicular./sum(n_perpendicular.^2)^0.5;
        n_parallel=n2(i,:)-n2(j,:); n_parallel=n_parallel./sum(n_parallel.^2)^0.5;
        
        r12=r(tip2V{Labels(i)}(1),:)-r(tip2V{Labels(j)}(1),:);
        %r_outofplane=abs(sum(r12.*n_outofplane));
        %r_perpendicular=abs(sum(r12.*n_perpendicular));
        %r_parallel=sum(r12.*n_parallel);
        D(i,j,4)=sum(r12.^2).^0.5;
        C(i,j,4)=sum(n2(i,:).*n2(j,:));
        O(i,j,4)=sum(r12.*n_parallel);
    end
end


% for i=1:length(Labels)
%     [e1,e2]=find(AMlbl==Labels(i));
%     e12=unique([e1;e2]);
%     AM_tree=AMlbl(e12,e12);
%     AM_tree(AM_tree~=Labels(i))=0;
%     AM_tree=((AM_tree+AM_tree')>0);
%     r_tree=r(e12,:);
%     
%     % disconnect tree branches
%     bp=find(sum(AM_tree)>2);
%     for j=1:length(bp)
%         bp_neigh=find(AM_tree(:,bp(j)));
%         r_tree=[r_tree;ones(length(bp_neigh)-1,1)*r_tree(bp(j),:)];
%         a=length(AM_tree);
%         AM_tree(a+length(bp_neigh)-1,a+length(bp_neigh)-1)=0;
%         AM_tree(bp_neigh(2:end),bp(j))=0;
%         AM_tree(bp(j),bp_neigh(2:end))=0;
%         temp_ind=sub2ind(size(AM_tree),bp_neigh(2:end),a-1+[2:length(bp_neigh)]');
%         AM_tree(temp_ind)=1;
%         temp_ind=sub2ind(size(AM_tree),a-1+[2:length(bp_neigh)]',bp_neigh(2:end));
%         AM_tree(temp_ind)=1;
%     end
% end





