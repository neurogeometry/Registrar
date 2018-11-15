% This function imports tif images into Matlab. RGB images are converted to
% grayscale. Stack can be imported in full or sparse data_format.
% Stacks can also be reduced in xy or xyz (reduct_type) by an integer
% factor, reduct_factor.
% Max, Min, Mean, Median, and SD reduct_method are implemented
% For the Full format, Orig is returned in uint8 or uint16 original data
% formats.
% For the Sparse format, Orig is double
% data_format is a string: 'Full' or 'Sparse'
% relative_thr is a double between 0 and 1
% reduct_type is a string: 'NA', 'xy', or 'xyz'
% reduct_factor is a double, e.g. 2
% reduct_method is a string: 'Min', 'Max', 'Mean', 'Median', or 'SD'

function [Orig,sizeOrig,classOrig]=ImportStack(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method)

b=[-0.0278,0.0319,0.0019];

if reduct_factor<1
    reduct_factor=1;
end
reduct_factor=fix(reduct_factor);

temp=dir([pth,'\*.tif']);
Names={temp.name};
for i=1:length(Names)
    Names{i}=Names{i}(1:find(Names{i}=='.')-1);
end
[~,ind]=sort(str2double(Names));
Names=Names(ind);
N=length(Names);

if strcmp(reduct_type,'NA')
    reduct_factor=1;
    reduct_amount=[1,1,1];
elseif strcmp(reduct_type,'xy')
    reduct_amount=[reduct_factor,reduct_factor,1];
elseif strcmp(reduct_type,'xyz')
    reduct_amount=reduct_factor.*[1,1,1];
end

temp = imread([pth,Names{i},'.tif']);
classOrig=class(temp);
formatOrig=size(temp,3);
sizeOrig=fix([size(temp,1),size(temp,2),N]./reduct_amount);
[xx,yy,zz]=ind2sub(reduct_amount,1:prod(reduct_amount));
if strcmp(data_format,'Full')
    Orig=zeros(sizeOrig,classOrig);
elseif strcmp(data_format,'Sparse')
    Orig = sparse(sizeOrig(1)*sizeOrig(2),sizeOrig(3));
end

if reduct_factor==1
    if strcmp(data_format,'Full')
        for i=1:N
            temp = imread([pth,Names{i},'.tif']);
            if formatOrig==3
                %temp = rgb2gray(temp);
                temp=temp(:,:,2);
                temp= intmax(classOrig)-temp;
            end
            Orig(:,:,i) = temp;
            display([num2str(i),'/',num2str(N)])
        end
    elseif strcmp(data_format,'Sparse')
        for i=1:N
            temp = imread([pth,Names{i},'.tif']);
            if formatOrig==3
                %temp = rgb2gray(temp);
                temp=temp(:,:,2);
                temp= intmax(classOrig)-temp;
            end
            %temp(temp<relative_thr*double(max(temp(:))))=0;
            temp(temp<relative_thr*double(intmax(classOrig)))=0;
            Orig(:,i)=sparse(double(temp(:)));
            display([num2str(i),'/',num2str(N)])
        end
        Orig=Orig(:);
    end
elseif reduct_factor>1
    for i=1:sizeOrig(3)
        Orig_temp=zeros([sizeOrig(1:2).*reduct_amount(1:2),reduct_amount(3)],classOrig);
        Orig3=zeros([sizeOrig(1:2),prod(reduct_amount)],classOrig);
        for j=1:reduct_amount(3)
            temp = imread([pth,Names{(i-1)*reduct_amount(3)+j},'.tif']);
            if formatOrig==3
                temp=intmax(classOrig)-temp;
                %%%%temp = rgb2gray(temp);
%                 temp=temp(:,:,2);
                                
                temp=double(temp);
                temp=temp(:,:,1)*b(1)+temp(:,:,2)*b(2)+temp(:,:,3)*b(3);
                temp=uint8((temp-min(temp(:)))./(max(temp(:))-min(temp(:))).*255);
            end
            Orig_temp(:,:,j)=temp(1:(sizeOrig(1)*reduct_amount(1)),1:(sizeOrig(2)*reduct_amount(2)));
        end
        for ii=1:prod(reduct_amount)
            Orig3(:,:,ii)=Orig_temp(xx(ii):reduct_amount(1):size(Orig_temp,1),yy(ii):reduct_amount(2):size(Orig_temp,2),zz(ii));
        end
        clear Orig_temp
        
        if strcmp(reduct_method,'Max')
            Orig3 = max(Orig3,[],3);
        elseif strcmp(reduct_method,'Min')
            Orig3 = min(Orig3,[],3);
        elseif strcmp(reduct_method,'Mean')
            Orig3 = mean(double(Orig3),3);
        elseif strcmp(reduct_method,'Median')
            Orig3 = median(double(Orig3),3);
        elseif strcmp(reduct_method,'SD')
            Orig3 = std(double(Orig3),[],3);
        else
            error('Unknown method')
        end
        if strcmp(data_format,'Full')
            Orig(:,:,i)=Orig3;
        elseif strcmp(data_format,'Sparse')
            Orig3=double(Orig3);
            Orig3(Orig3<relative_thr*double(intmax(classOrig)))=0;
            Orig(:,i)=sparse(Orig3(:));
        end
        display([num2str(i),'/',num2str(sizeOrig(3))])
    end
    if strcmp(data_format,'Sparse')
        Orig=Orig(:);
    end
end
disp('Stack is imported.')