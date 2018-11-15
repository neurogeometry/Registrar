% This function performs Max, Min, Mean, Median, and SD filtering
% filter_type is a string: 'Min', 'Max', 'Mean', 'Median', or 'SD'
% filter_size is a vector of three odd integers, e.g. [3,3,1]

function Im1=Simple_Filters(Im,filter_type,filter_size)

filter_size=fix(filter_size./2)*2+1;
disp([num2str(filter_size),' ',filter_type,' filtering started.'])

cl=class(Im);
sizeIm=size(Im);
if length(sizeIm)==2
    sizeIm(3)=1;
end
Im=double(Im);
Im1=zeros(sizeIm);
HW=(filter_size-1)./2;

for i=1:sizeIm(3)
    beg_ind=max(1,i-HW(3));
    end_ind=min(sizeIm(3),i+HW(3));
    Im_temp=Im(:,:,beg_ind:end_ind);
        
    if strcmp(filter_type,'Max')
        Im1(:,:,i) = max(Im_temp,[],3);
        Im1(:,:,i) = ordfilt2(Im1(:,:,i), filter_size(1)*filter_size(2), ones(filter_size(1),filter_size(2)));
    elseif strcmp(filter_type,'Min')
        Im1(:,:,i) = min(Im_temp,[],3);
        Im1(:,:,i) = ordfilt2(Im1(:,:,i), 1, ones(filter_size(1),filter_size(2)));
    elseif strcmp(filter_type,'Mean')
        Im1(:,:,i) = mean(Im_temp,3);
        H=fspecial('average',filter_size(1:2));
        Im1(:,:,i) = imfilter(Im1(:,:,i),H,'replicate');
    elseif strcmp(filter_type,'Median') % this calculation does not give the true median
        Im1(:,:,i) = median(Im_temp,3); 
        Im1(:,:,i) = ordfilt2(Im1(:,:,i), (filter_size(1)*filter_size(2)+1)/2, ones(filter_size(1),filter_size(2)));
    elseif strcmp(filter_type,'SD')
        X2 = mean(Im_temp.^2,3);
        X1 = mean(Im_temp,3);
        H=fspecial('average',filter_size(1:2));
        Im1(:,:,i) = (imfilter(X2,H,'replicate')-imfilter(X1,H,'replicate').^2).^0.5;
    else
        error('Unknown method')
    end 
end

if ~strcmp(cl,'double')
    Im1=Im1-min(Im1(:));
    Im1=feval(cl,Im1./max(Im1(:)).*double(intmax(cl)));
end
disp([num2str(filter_size),' ',filter_type,' filtering is complete'])