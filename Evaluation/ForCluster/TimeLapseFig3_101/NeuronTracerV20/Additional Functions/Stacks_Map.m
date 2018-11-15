% This function imports and displays multiple stacks of images by using ImportStack.m
% Approximate initial x,y positions of all the stacks must be provided in 
% Initial_Stack_xy_Positions.txt file located in the same folder.
% Initial_Start_xy_Points.txt or Initial_Start_xyz_Points.txt contains the start point pasitions (if available)
% ref_stack is the reference frame for the start points

function Stacks_Map(pth,ref_stack)

data_format='Full'; 
relative_thr=0.1;
reduct_type='NA'; % 'NA','xy','xyz'
reduct_factor=1;
reduct_method='Median';

if reduct_factor<1
    reduct_factor=1;
end
reduct_factor=fix(reduct_factor);

if exist([pth,'Initial_Start_xy_Points.txt'],'file')
    SP_xy=textread([pth,'Initial_Start_xy_Points.txt']); % java format
    SP_xy=fix(SP_xy(:,[2,1])./reduct_factor)+1; % matlab format
elseif exist([pth,'Initial_Start_xyz_Points.txt'],'file')
    SP_xyz=textread([pth,'Initial_Start_xyz_Points.txt']);
    SP_xy=fix(SP_xyz(:,[2,1])./reduct_factor)+1; 
end

if exist([pth,'Initial_Stack_xyz_Positions.txt'],'file')
    Stack_xy_positions=textread([pth,'Initial_Stack_xyz_Positions.txt']); % java format
    Stack_xy_positions=fix(Stack_xy_positions(:,[2,1])./reduct_factor)+1; % matlab format
elseif exist([pth,'Initial_Stack_xy_Positions.txt'],'file')
    Stack_xy_positions=textread([pth,'Initial_Stack_xy_Positions.txt']); % java format
    Stack_xy_positions=fix(Stack_xy_positions(:,[2,1])./reduct_factor)+1; % matlab format
end


temp=dir(pth);
StackNames={temp.name};
StackNames=StackNames([temp.isdir]);
StackNames=StackNames(3:end);
StackNames(strcmp(StackNames,'Results'))=[]; %exclude the Results folder
[temp,ind]=sort(str2double(StackNames));
StackNames=StackNames(ind);
N_stacks=length(StackNames);

% import stacks
sizeStack=zeros(N_stacks,3);
for i=1:N_stacks,i
    [temp,sizeStack(i,:),classStack]=ImportStack([pth,StackNames{i},'\'],data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
    temp=255-temp;
    Stack{i}=max(temp,[],3);
end

sizeIm=[max(Stack_xy_positions+sizeStack(:,[1,2]))-min(Stack_xy_positions)];
Im=feval(classStack,zeros(sizeIm));
for i=1:N_stacks
    temp1=Stack_xy_positions(i,:)-min(Stack_xy_positions)+1;
    
    Im(temp1(1):temp1(1)+sizeStack(i,1)-1,temp1(2):temp1(2)+sizeStack(i,2)-1)=max(Im(temp1(1):temp1(1)+sizeStack(i,1)-1,temp1(2):temp1(2)+sizeStack(i,2)-1),Stack{i});
end

figure 
imshow(Im,[0 max(Im(:))])
%imshow(Im,[0 2^12-1])
hold on
title('Mosaic')
for i=1:N_stacks
    text(Stack_xy_positions(i,2)-min(Stack_xy_positions(:,2))+1,Stack_xy_positions(i,1)-min(Stack_xy_positions(:,1))+1,num2str(i),'Color','r')
end

if isempty(ref_stack)
    for i=1:size(SP_xy,1)
        plot(SP_xy(i,2)-min(Stack_xy_positions(:,2)),SP_xy(i,1)-min(Stack_xy_positions(:,1)),'c*')
        text(SP_xy(i,2)-min(Stack_xy_positions(:,2))+5,SP_xy(i,1)-min(Stack_xy_positions(:,1))+5,num2str(i),'Color','c')
    end
else
    for i=1:size(SP_xy,1)
        plot(SP_xy(i,2)+Stack_xy_positions(ref_stack,2)-1,SP_xy(i,1)+Stack_xy_positions(ref_stack,1)-1,'c*')
        text(SP_xy(i,2)+Stack_xy_positions(ref_stack,2)-1+5,SP_xy(i,1)+Stack_xy_positions(ref_stack,1)-1+5,num2str(i),'Color','r')
    end
end

%imwrite(Im,[pth,'Results\map.tif'],'tif')
