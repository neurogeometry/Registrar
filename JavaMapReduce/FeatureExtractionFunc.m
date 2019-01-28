function FeatureExtractionFunc(File,stackID,DataFolder,StackPositions_pixels,StackSizes_pixels)
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
% Northeastern University, USA
% =========================================================================
% -------------------------------------------------------------------------
FEwindowsize = [4 4 4];

% Import Stack
IM_Original=ImportStack(char(File),StackSizes_pixels(1,:));

% Find seeds
r_seed=Find_Seeds(IM_Original,StackPositions_pixels,StackSizes_pixels);

% Remove extra features
[X1,Y1,Z1]=size(IM_Original);
% if Z1 < FEwindowsize(3)
%     FEwindowsize(3) = 0;
% end
rem_ind_1 = (r_seed(:,1) <=FEwindowsize(1) | r_seed(:,1) > (X1-FEwindowsize(1)) | r_seed(:,2) <=FEwindowsize(2) | r_seed(:,2) > (Y1-FEwindowsize(2)) | r_seed(:,3) <=FEwindowsize(3) | r_seed(:,3) > (Z1-FEwindowsize(3)));
r_seed(rem_ind_1,:) = [];

% Get Feature vectors
FeatureVector=[];
tic
FeatureVector = zeros(size(r_seed,1),prod(2*FEwindowsize+1));
for i=1:size(r_seed,1)
    temp = double(IM_Original(r_seed(i,1)-FEwindowsize(1):r_seed(i,1)+FEwindowsize(1),r_seed(i,2)-FEwindowsize(2):r_seed(i,2)+FEwindowsize(2),r_seed(i,3)-FEwindowsize(3):r_seed(i,3)+FEwindowsize(3)));
    FeatureVector(i,:) =  temp(:)';
end

% Write feature locations and Feature vectors into HDF5 files
hdf5write([DataFolder,'/tmp/Feature_seeds',num2str(stackID),'.h5'], '/dataset1', r_seed);
hdf5write([DataFolder,'/tmp/Feature_vector',num2str(stackID),'.h5'], '/dataset1', FeatureVector);
seedsFile = [DataFolder,'/tmp/Feature_',num2str(stackID),'.h5'];
end

