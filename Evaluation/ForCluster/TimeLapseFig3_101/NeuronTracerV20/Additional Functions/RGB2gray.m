% This function converts GRB to grayscale format in a way that optimizes a
% particular cost function.

function Orig=RGB2gray(pth)

nn=10; % sample image
Nblobs=20;
Naxons=20;
Nbackg=20;

temp=dir([pth,'\*.tif']);
Names={temp.name};
for i=1:length(Names)
    Names{i}=Names{i}(1:find(Names{i}=='.')-1);
end
[~,ind]=sort(str2double(Names));
Names=Names(ind);
N=length(Names);

temp=imread([pth,Names{nn},'.tif']);
classOrig=class(temp);
formatOrig=size(temp,3);
sizeOrig=fix([size(temp,1),size(temp,2),N]);

temp=(255-temp);
temp_gray=double(rgb2gray(temp));
temp=double(temp);
R=temp(:,:,1);
G=temp(:,:,2);
B=temp(:,:,3);

figure, hold on
imshow(temp_gray,[0 max(temp_gray(:))])
[x,y]=ginput(Nblobs);
ind=sub2ind(size(temp_gray),fix(y),fix(x));
RGB_blobs=[R(ind),G(ind),B(ind),ones(length(ind),1)];
[x,y]=ginput(Naxons);
ind=sub2ind(size(temp_gray),fix(y),fix(x));
RGB_axons=[R(ind),G(ind),B(ind),ones(length(ind),1)];
[x,y]=ginput(Nbackg);
ind=sub2ind(size(temp_gray),fix(y),fix(x));
RGB_backg=[R(ind),G(ind),B(ind),ones(length(ind),1)];

% R_relative=RGB_axons(:,1)*ones(1,Nblobs)-ones(Naxons,1)*RGB_blobs(:,1)';
% G_relative=RGB_axons(:,2)*ones(1,Nblobs)-ones(Naxons,1)*RGB_blobs(:,2)';
% B_relative=RGB_axons(:,3)*ones(1,Nblobs)-ones(Naxons,1)*RGB_blobs(:,3)';
% RGB_relative=[R_relative(:),G_relative(:),B_relative(:)];

b=(RGB_axons'*RGB_axons./Naxons+1.*RGB_blobs'*RGB_blobs./Nblobs+1.*RGB_backg'*RGB_backg./Nbackg)^(-1)*(mean(RGB_axons).*1+mean(RGB_blobs).*(0)+mean(RGB_backg).*(0))';

% [B,lamda]=eig(RGB_relative'*RGB_relative);
% b=B(:,3);

Orig=zeros(sizeOrig,classOrig);
for i=1:N, i/N
    temp = 255-double(imread([pth,Names{i},'.tif']));
    if formatOrig==3
        temp=temp(:,:,1)*b(1)+temp(:,:,2)*b(2)+temp(:,:,3)*b(3)+ones(size(temp(:,:,1)))*b(4);
        temp=uint8((temp-min(temp(:)))./(max(temp(:))-min(temp(:))).*255);
    end
    Orig(:,:,i) = temp;
end
%Orig=uint8((Orig-min(Orig(:)))./(max(Orig(:))-min(Orig(:))).*255);
