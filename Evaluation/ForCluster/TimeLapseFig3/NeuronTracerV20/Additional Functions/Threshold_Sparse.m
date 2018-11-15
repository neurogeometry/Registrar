% This function thresholds the stack in sparse format.

function Im=Threshold_Sparse(Im,thr_low)

[ind,~,v]=find(Im);
v(v<thr_low)=0;
Im=sparse(ind,1,v,length(Im),1);
