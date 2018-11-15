% This function disconnects all branches at branch points
% AMlbl in the input must be labeled with LabelBranchesAM
% step back is in the amount of min(segment_length/2, 1)

function [AMlbl,r,R] = Disconnect_Branches(AMlbl,r,R)

del=0.5;

AMlbl=LabelBranchesAM(AMlbl);
bps=find(sum((AMlbl>0))>=3);
for i=1:length(bps)
    bp_labels=nonzeros(AMlbl(:,bps(i)));
    if length(bp_labels)>1
        r=[r;ones(length(bp_labels)-1,1)*r(bps(i),:)];
        R=[R;ones(length(bp_labels)-1,1)*R(bps(i))];
        a=length(AMlbl);
        AMlbl(a+length(bp_labels)-1,a+length(bp_labels)-1)=0;
        
        temp=find(AMlbl(:,bps(i))==bp_labels(1),1,'first');
        del_temp=min(del,1/sum((r(temp,:)-r(bps(i),:)).^2)^0.5);
        r(bps(i),:)=r(bps(i),:)+(r(temp,:)-r(bps(i),:)).*del_temp;
        R(bps(i))=R(bps(i))+(R(temp)-R(bps(i))).*del_temp;
        for j=2:length(bp_labels)
            temp=find(AMlbl(:,bps(i))==bp_labels(j),1,'first');
            AMlbl(temp,bps(i))=0;
            AMlbl(bps(i),temp)=0;
            AMlbl(temp,a+j-1)=bp_labels(j);
            AMlbl(a+j-1,temp)=bp_labels(j);
            del_temp=min(del,1/sum((r(temp,:)-r(a+j-1,:)).^2)^0.5);
            r(a+j-1,:)=r(a+j-1,:)+(r(temp,:)-r(a+j-1,:)).*del_temp;
            R(a+j-1)=R(a+j-1)+(R(temp)-R(a+j-1)).*del_temp;
        end
    end
end
%AMlbl = LabelBranchesAM(spones(AMtemp));
