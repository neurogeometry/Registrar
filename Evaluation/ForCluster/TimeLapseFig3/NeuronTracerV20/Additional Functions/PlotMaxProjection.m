% This function shows the stack max projection, from file or Im.

%not finished !!!!!!!!!!!!!
function PlotMaxProjection(pth,Im,sizeIm,data_format,relative_thr,reduct_type,reduct_factor,reduct_method)

if ~isempty(pth)
    [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method);

elseif ~isempty(Im)
    MaxOrig=zeros(sizeIm(1:2));
    [ind,~,v]=find(Im);
    if ~isempty(ind)
        [i,j,k]=ind2sub_AS(sizeIm,ind);
        ind_ij=sub2ind(sizeIm(1:2),i,j);
        [ind_ij,sort_ind]=sort(ind_ij);
        v=v(sort_ind);
        temp=[find((ind_ij(2:end)-ind_ij(1:end-1))>0);length(ind_ij)];
        ind_ij=ind_ij(temp);
        temp_max=zeros(length(temp),1);
        temp_max(1)=max(v(1:temp(1)));
        for ii=2:length(temp)
            temp_max(ii)=max(v(temp(ii-1):temp(ii)));
        end
        MaxOrig(ind_ij)=temp_max;
    end
end

imshow(MaxOrig,[0 max(MaxOrig(:))]);

