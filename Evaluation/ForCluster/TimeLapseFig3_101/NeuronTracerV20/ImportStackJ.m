% This function imports images into MatLab. RGB images are converted to grayscale.
% Stacks can be reduced in xy and z by integer factors (reductions).
% Max, Min, Mean, Median, and SD reduct_method are implemented.
% Data format of Orig is preserved (uint8, uint16, etc).
% reduction_x, reduction_y, and reduction_z are doubles, e.g. 2.1
% reduct_method is a string: 'Min', 'Max', 'Mean', 'Median', or 'SD'

function [Orig,sizeOrig,classOrig]=ImportStackJ(pth,file_list,reduction_x,reduction_y,reduction_z,reduct_method)

Orig=[];
sizeOrig=[];
classOrig=[];

N=length(file_list);
info = imfinfo([pth,file_list{1}]);
Npl=length(info);

reduction_y=round(reduction_y);
if reduction_y<1
    reduction_y=1;
end
reduction_x=round(reduction_x);
if reduction_x<1
    reduction_x=1;
end
reduction_z=round(reduction_z);
if reduction_z<1
    reduction_z=1;
end
    
if N==0
    disp('There are no images which pass the file selection criteria.')
    return  
elseif N==1 && Npl>1 % import a virtual stack (tif or LSM)
    
    temp = imread([pth,file_list{1}],'Index',1);
    
    if reduction_z>Npl
        reduction_z=Npl;
    end
    reduct_amount=[reduction_y,reduction_x,reduction_z];
    
    classOrig=class(temp);
    formatOrig=size(temp,3);
    sizeOrig=fix([size(temp,1),size(temp,2),Npl]./reduct_amount);
    [xx,yy,zz]=ind2sub(reduct_amount,1:prod(reduct_amount));
    Orig=zeros(sizeOrig,classOrig);
    
    if reduction_y==1 && reduction_x==1 && reduction_z==1
        for i=1:Npl
            temp = imread([pth,file_list{1}],'Index',i);
            if formatOrig==3
                temp=rgb2gray(temp);
            end
            Orig(:,:,i) = temp;
            display([num2str(i),'/',num2str(Npl)])
        end
    elseif reduction_y>1 || reduction_x>1 || reduction_z>1
        for i=1:sizeOrig(3)
            Orig_temp=zeros([sizeOrig(1:2).*reduct_amount(1:2),reduct_amount(3)],classOrig);
            Orig3=zeros([sizeOrig(1:2),prod(reduct_amount)],classOrig);
            for j=1:reduct_amount(3)
                temp = imread([pth,file_list{1}],'Index',(i-1)*reduct_amount(3)+j);
                if formatOrig==3
                    temp=rgb2gray(temp);
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
                Orig3=sort(Orig3,3);
                if mod(size(Orig3,3),2)==0
                    Orig3=mean(Orig3(:,:,size(Orig3,3)/2:size(Orig3,3)/2+1),3);
                else
                    Orig3=Orig3(:,:,(size(Orig3,3)+1)/2);
                end
                %Orig3 = median(double(Orig3),3);
            elseif strcmp(reduct_method,'SD')
                Orig3 = std(double(Orig3),[],3);
            else
                error('Unknown method')
            end
            Orig(:,:,i)=Orig3;
            display([num2str(i),'/',num2str(sizeOrig(3))])
        end
    end    
elseif N>1 || (N==1 && Npl==1) % import a set of regular images
    
    if Npl>1
        disp('Unable to load multiple virtual stacks.')
        return
    end
    
    temp = imread([pth,file_list{1}]);
    for i=2:N
        info = imfinfo([pth,file_list{i}]);
        if length(info)>1
            disp('Unable to load multiple virtual stacks.')
            return
        end
        if info.Height~=size(temp,1) || info.Width~=size(temp,2)
            disp('Unable to load. Images have different dimensions.')
            return
        end
    end
    
    if reduction_z>N
        reduction_z=N;
    end
    reduct_amount=[reduction_y,reduction_x,reduction_z];
    
    classOrig=class(temp);
    formatOrig=size(temp,3);
    sizeOrig=fix([size(temp,1),size(temp,2),N]./reduct_amount);
    [xx,yy,zz]=ind2sub(reduct_amount,1:prod(reduct_amount));
    Orig=zeros(sizeOrig,classOrig);
    
    if reduction_y==1 && reduction_x==1 && reduction_z==1
        for i=1:N
            temp = imread([pth,file_list{i}]);
            if formatOrig==3
                temp=rgb2gray(temp);
            end
            Orig(:,:,i) = temp;
            display([num2str(i),'/',num2str(N)])
        end
    elseif reduction_y>1 || reduction_x>1 || reduction_z>1
        for i=1:sizeOrig(3)
            Orig_temp=zeros([sizeOrig(1:2).*reduct_amount(1:2),reduct_amount(3)],classOrig);
            Orig3=zeros([sizeOrig(1:2),prod(reduct_amount)],classOrig);
            for j=1:reduct_amount(3)
                temp = imread([pth,file_list{(i-1)*reduct_amount(3)+j}]);
                if formatOrig==3
                    temp=rgb2gray(temp);
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
                Orig3=sort(Orig3,3);
                if mod(size(Orig3,3),2)==0
                    Orig3=mean(Orig3(:,:,size(Orig3,3)/2:size(Orig3,3)/2+1),3);
                else
                    Orig3=Orig3(:,:,(size(Orig3,3)+1)/2);
                end
                %Orig3 = median(double(Orig3),3);
            elseif strcmp(reduct_method,'SD')
                Orig3 = std(double(Orig3),[],3);
            else
                error('Unknown method')
            end
            Orig(:,:,i)=Orig3;
            display([num2str(i),'/',num2str(sizeOrig(3))])
        end
    end
    disp('Stack is imported.')
end

