% This function shifts and combines SWC files of multiple stacks
% Approximate initial x,y positions of all the stacks must be provided in 
% Initial_Stack_xy_Positions.txt file located in the same folder.

function Combine_SWC(pth,Stack_size)

if exist([pth,'Initial_Stack_xyz_Positions.txt'],'file')
    Stack_xyz_positions=textread([pth,'Initial_Stack_xyz_Positions.txt']); % java format
elseif exist([pth,'Initial_Stack_xy_Positions.txt'],'file')
    Stack_xyz_positions=textread([pth,'Initial_Stack_xy_Positions.txt']); % java format
    Stack_xyz_positions=[Stack_xyz_positions,zeros(size(Stack_xyz_positions,1),1)];
end

temp=dir(pth);
StackNames={temp.name};
StackNames=StackNames([temp.isdir]);
StackNames=StackNames(3:end);
StackNames(strcmp(StackNames,'Results'))=[]; %exclude the Results folder
[temp,ind]=sort(str2double(StackNames));
StackNames=StackNames(ind);
N_stacks=length(StackNames);

% import SWC files
SWC=[];
AMall=sparse([]);
Rall=[];
for i=1:N_stacks,i
%     load([pth,'Results\',StackNames{i},'_AM-R.mat']);
%     AM(AM~=0)=AM(AM~=0)+10^5*i;
%     R=R+ones(size(R,1),1)*(Stack_xyz_positions(i,[2,1,3])-min(Stack_xyz_positions(:,[2,1,3])));
%     
%     %trim
%     rem_ind=[];
%     for j=i+1:N_stacks
%         rem_temp=find(R(:,1)>Stack_xyz_positions(j,2)-min(Stack_xyz_positions(:,2)) & ...
%             R(:,1)<Stack_xyz_positions(j,2)-min(Stack_xyz_positions(:,2))+Stack_size & ...
%             R(:,2)>Stack_xyz_positions(j,1)-min(Stack_xyz_positions(:,1)) & ...
%             R(:,2)<Stack_xyz_positions(j,1)-min(Stack_xyz_positions(:,1))+Stack_size);
%         
%         rem_ind=[rem_ind;rem_temp];
%     end
%     AM(:,rem_ind)=[];
%     AM(rem_ind,:)=[];
%     R(rem_ind,:)=[];
%     AM=LabelBranchesAM_AS(AM~=0);
%     AM(AM~=0)=AM(AM~=0)+10^5*i;
%     
%     AMall(end+(1:size(R,1)),end+(1:size(R,1)))=AM;
%     Rall=[Rall;R];
    
    SWCtemp=textread([pth,'Results\',StackNames{i},'_All.swc']);
    SWCtemp(:,8)=[];
    SWCtemp(:,1)=SWCtemp(:,1)+size(SWC,1);
    SWCtemp(:,3)=SWCtemp(:,3)+Stack_xyz_positions(i,1);
    SWCtemp(:,4)=SWCtemp(:,4)+Stack_xyz_positions(i,2);
    SWCtemp(:,5)=SWCtemp(:,5)+Stack_xyz_positions(i,3);
    SWCtemp(SWCtemp(:,7)~=-1,7)=SWCtemp(SWCtemp(:,7)~=-1,7)+size(SWC,1);
    SWC=[SWC;SWCtemp];
end
%[~, SWC]=AM2swc(AMall,Rall,'NA',1);

dlmwrite([pth,'Results\','6_Section4_Combined.swc'], SWC, 'delimiter', ' ')



