% This function finds clusters in a directed or an undirected AM.
% ClusterV is a cell array of cluster vertices.
% ClusterN is the number of clusters.

function [ClusterN,ClusterV]=FindClustersAM(AM)

AM=logical((AM+AM'));
ClusterN=0;
ClusterV={};

[~,V1]=find(AM,1,'first');
if ~isempty(V1)
    ClusterN=1;
    V=false(1,size(AM,1));
    V(V1)=true;
    ClusterV{ClusterN}=V; 
end

while ~isempty(V1)
    Vnew=logical(sum(AM(find(V),:),1)); % fastest way for sparse AM
    if sum(Vnew)>0
        ClusterV{ClusterN}=(ClusterV{ClusterN} | Vnew);
        AM(V,Vnew)=0;
        AM(Vnew,V)=0;
        V=Vnew;
    else
        ClusterV{ClusterN}=find(ClusterV{ClusterN});
        [~,V1]=find(AM,1,'first');
        if ~isempty(V1)
            ClusterN=ClusterN+1;
            V=false(1,size(AM,1));
            V(V1)=true;
            ClusterV{ClusterN}=V;
        end
    end
end
