% This function finds branch lengths in a labeled undirected AMlbl

function BranchLengths=BranchLengthsAM(AMlbl,r)

AMlbl_direct=triu(AMlbl,1);
L=unique(AMlbl_direct(AMlbl_direct>0));
BranchLengths=zeros(1,max(L));

for i=1:length(L)
    [e1 e2]=find(AMlbl==L(i));
    BranchLengths(L(i))=sum(sum((r(e1,:)-r(e2,:)).^2,2).^0.5)/2;
end

%disp(BranchLengths)