% This function determines if there are loops in a directed or an undirected AM

function YN = Loops(AM)

AM=spones(AM+AM'); AM=AM-diag(diag(AM));
temp_ind=(sum(AM)~=0);
AM=AM(temp_ind,temp_ind);

AV=find(sum(AM));
NV=length(AV);
if NV==0
    YN=0;
else
    startV=AV(1);
    visited=startV;
    AMV=AM(startV,:);
    startVnew=find(sum(AMV,1));
    while ~isempty(AV) && nnz(AMV)==length(startVnew)
        if ~isempty(startVnew)
            AM(startV,startVnew)=0;
            AM(startVnew,startV)=0;
            startV=startVnew;
            visited=[visited,startV];
        else
            AV=find(sum(AM));
            if ~isempty(AV)
                startV=AV(1);
                visited=[visited,startV];
            end
        end
        AMV=AM(startV,:);
        startVnew=find(sum(AMV,1));
    end
    if length(visited)~=NV
        YN=1;
    else
        YN=0;
    end
end



