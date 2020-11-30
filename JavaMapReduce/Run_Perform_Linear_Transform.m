

FilePath = 'E:\Datasets\DIADEM\Neocortical Layer 1 Axons\Subset1_c\00001\00001.tif';
StackSize = [512,512,60];
StackPosition = [1,1,1];

L = [0.86,0.50,0;
    0.50,0.86,0;
    0,0,1];
b = [1;1;1];

X=ImportStack(FilePath,StackSize);
figure,
imshow(max(X,[],3),[0 max(X(:))]);
title('Original');
[X_prime,StackPosition_prime]=Perform_Linear_Transform(X,StackPosition,L,b);
figure,
imshow(max(X_prime,[],3),[0 max(X_prime(:))]);
title('Transformed');