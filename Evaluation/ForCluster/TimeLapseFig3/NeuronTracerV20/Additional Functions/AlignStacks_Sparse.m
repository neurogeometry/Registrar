% This function imports multiple stacks of images by using ImportStack.m
% and performs a rigid shift alignment. 
% Stacks are first aligned in xy, then in z, and finally in xyz. The result
% is saved in a sparse double formats as a single stack. 
% Approximate initial x,y positions of all the stacks must be provided in 
% Initial_Mosaic_Configuration.txt file located in the same folder. 

function [Orig,sizeOrig]=AlignStacks_Sparse(pth,data_format,relative_thr,reduct_type,reduct_factor,reduct_method)

addpath('C:\Armen\DIADEM\Neuron Tracer V2')
if reduct_factor<1
    reduct_factor=1;
end
reduct_factor=fix(reduct_factor);

if exist([pth,'Initial_Start_xy_Point.txt'],'file')
    SP_xy=textread([pth,'Initial_Start_xy_Point.txt']);
    SP_xy=fix(SP_xy(:,[2,1])./reduct_factor)+1; %%% !!!!!
elseif exist([pth,'Initial_Start_xyz_Point.txt'],'file')
    SP_xyz=textread([pth,'Initial_Start_xyz_Point.txt']);
    SP_xy=fix(SP_xyz(:,[2,1])./reduct_factor)+1; %%% !!!!!
end
Stack_xy_positions=textread([pth,'Initial_Stack_xy_Positions.txt']);
Stack_xy_positions=fix(Stack_xy_positions(:,[2,1])./reduct_factor)+1; %%% !!!!!

temp=dir(pth);
StackNames={temp.name};
StackNames=StackNames([temp.isdir]);
StackNames=StackNames(3:end);
[temp,ind]=sort(str2double(StackNames));
StackNames=StackNames(ind);
N_stacks=length(StackNames);

% import stacks
sizeStack=zeros(N_stacks,3);
for i=1:N_stacks,i
    [Stack{i},sizeStack(i,:)]=ImportStack([pth,StackNames{i},'\'],data_format,relative_thr,reduct_type,reduct_factor,reduct_method);
    whos Stack
end

% for Neuromuscular stacks it is better to save at this point
% save([pth,'Original_Mosaic_3Median_0_01.mat'],'Stack','sizeStack')
% load([pth,'Original_Mosaic_3Median_0_005.mat']);
% Stack{91}=Stack{91}./2^4; % !!!
% Stack{92}=Stack{92}./2^4; % !!!


%Stack=Stack(100:end); Stack_xy_positions=Stack_xy_positions(100:end,:); sizeStack=sizeStack(100:end,:); N_stacks=length(Stack);

% combine stacks and desplay max projection
[Orig,sizeOrig]=Combine_Stacks(Stack,sizeStack,[Stack_xy_positions,ones(N_stacks,1)]);
figure 
PlotMaxProjection_Sparse([],Orig,sizeOrig)
hold on
title('Original Mosaic')
for i=1:N_stacks
    text(Stack_xy_positions(i,2)-min(Stack_xy_positions(:,2))+1,Stack_xy_positions(i,1)-min(Stack_xy_positions(:,1))+1,num2str(i),'Color','r')
end
plot(SP_xy(:,2)-min(Stack_xy_positions(:,2))+1,SP_xy(:,1)-min(Stack_xy_positions(:,1))+1,'r*')

% align stack projections in xy
for i=1:N_stacks,i
    StackProj{i}=MaxProjection(Stack{i},sizeStack(i,:));
end

thr=50;
AlignRegion=[-15:15];
Cmax=zeros(N_stacks,N_stacks);
indmax=zeros(N_stacks,N_stacks);
Nmax=zeros(N_stacks,N_stacks);
[Align_x,Align_y]=meshgrid(AlignRegion,AlignRegion);
for i=1:N_stacks,i
    for j=(i+1):N_stacks
        if Stack_xy_positions(i,1)+sizeStack(i,1)>=Stack_xy_positions(j,1) && Stack_xy_positions(j,1)+sizeStack(j,1)>=Stack_xy_positions(i,1) && Stack_xy_positions(i,2)+sizeStack(i,2)>=Stack_xy_positions(j,2) && Stack_xy_positions(j,2)+sizeStack(j,2)>=Stack_xy_positions(i,2)
            Ctemp=zeros(1,numel(Align_x));
            Ntemp=zeros(size(Ctemp));
            nn=0;
            for k=1:numel(Align_x)
                nn=nn+1;
                [Ctemp(nn),Ntemp(nn)] = Corr_Sparse_xy(StackProj{i},StackProj{j},sizeStack(i,1:2),sizeStack(j,1:2),Stack_xy_positions(i,:)+[Align_x(k),Align_y(k)],Stack_xy_positions(j,:),thr);
            end
            [~,indmax(i,j)]=max(Ctemp(:).*Ntemp(:).^0);
            Cmax(i,j)=Ctemp(indmax(i,j));
            Nmax(i,j)=Ntemp(indmax(i,j));
        end
    end
end
x_max=zeros(size(indmax));
x_max(indmax>0)=Align_x(indmax(indmax>0));
y_max=zeros(size(indmax));
y_max(indmax>0)=Align_y(indmax(indmax>0));

Cmax=Cmax+Cmax';
Nmax=Nmax+Nmax';
x_max=x_max-x_max';
y_max=y_max-y_max';
W=Cmax.*Nmax.^0.5; % weight

Stack_xy_aligned_positions = Stack_xy_positions+round((W-diag(sum(W)-10^-6))^(-1)*[diag(W*x_max),diag(W*y_max)]);
Stack_xy_aligned_positions=[Stack_xy_aligned_positions(:,1)-min(Stack_xy_aligned_positions(:,1)),Stack_xy_aligned_positions(:,2)-min(Stack_xy_aligned_positions(:,2))]+1;


% Stack_xy_aligned_positions=Stack_xy_positions;
% for k=1:50
%     C_xy=zeros(N_stacks,N_stacks);
%     Ntemp=zeros(N_stacks,N_stacks);
%     C_xp1y=zeros(N_stacks,N_stacks);
%     C_xm1y=zeros(N_stacks,N_stacks);
%     C_xyp1=zeros(N_stacks,N_stacks);
%     C_xym1=zeros(N_stacks,N_stacks);
%     for i=1:N_stacks
%         for j=(i+1):N_stacks
%             [C_xy(i,j),Ntemp(i,j)] = Corr_Sparse_xy(StackProj{i},StackProj{j},sizeStack(i,1:2),sizeStack(j,1:2),Stack_xy_aligned_positions(i,:),Stack_xy_aligned_positions(j,:),thr);
%             [C_xp1y(i,j),~] = Corr_Sparse_xy(StackProj{i},StackProj{j},sizeStack(i,1:2),sizeStack(j,1:2),Stack_xy_aligned_positions(i,:)+[1,0],Stack_xy_aligned_positions(j,:),thr);
%             [C_xm1y(i,j),~] = Corr_Sparse_xy(StackProj{i},StackProj{j},sizeStack(i,1:2),sizeStack(j,1:2),Stack_xy_aligned_positions(i,:)+[-1,0],Stack_xy_aligned_positions(j,:),thr);
%             [C_xyp1(i,j),~] = Corr_Sparse_xy(StackProj{i},StackProj{j},sizeStack(i,1:2),sizeStack(j,1:2),Stack_xy_aligned_positions(i,:)+[0,1],Stack_xy_aligned_positions(j,:),thr);
%             [C_xym1(i,j),~] = Corr_Sparse_xy(StackProj{i},StackProj{j},sizeStack(i,1:2),sizeStack(j,1:2),Stack_xy_aligned_positions(i,:)+[0,-1],Stack_xy_aligned_positions(j,:),thr);
%         end
%     end
%     W=(Ntemp.^0.5); W=W+W';
%     %W=1;
%     tempx=sum(((C_xp1y+C_xm1y')-(C_xm1y+C_xp1y')).*W,2);
%     %[tempxmax,imax]=max(abs(tempx));
%     tempy=sum(((C_xyp1+C_xym1')-(C_xym1+C_xyp1')).*W,2);
%     %[tempymax,jmax]=max(abs(tempy));
%     %Stack_xy_aligned_positions(imax,1)=Stack_xy_aligned_positions(imax,1)+sign(tempxmax);
%     %Stack_xy_aligned_positions(jmax,2)=Stack_xy_aligned_positions(jmax,2)+sign(tempymax);
%     tempx(abs(tempx)<0.2.*max(abs([tempx;tempy])))=0;
%     tempy(abs(tempy)<0.2.*max(abs([tempx;tempy])))=0;
%     Stack_xy_aligned_positions=Stack_xy_aligned_positions+sign([tempx,tempy]);
%     [k,mean(C_xy(:).*W(:))]
% end

[Orig,sizeOrig]=Combine_Stacks(Stack,sizeStack,[Stack_xy_aligned_positions,ones(length(Stack),1)]);
figure
PlotMaxProjection_Sparse([],Orig,sizeOrig)
title('Mosaic after xy alignment')
for i=1:N_stacks
    text(Stack_xy_aligned_positions(i,2)-min(Stack_xy_aligned_positions(:,2))+1,Stack_xy_aligned_positions(i,1)-min(Stack_xy_aligned_positions(:,1))+1,num2str(i),'Color','r')
end

% align stacks in z
Stack_z_aligned_positions=ones(N_stacks,1);

Cmax=zeros(N_stacks,N_stacks);
Nmax=zeros(N_stacks,N_stacks);
delZmax=zeros(N_stacks,N_stacks);
for i=1:N_stacks,i
    for j=(i+1):N_stacks
        if Stack_xy_aligned_positions(i,1)+sizeStack(i,1)>=Stack_xy_aligned_positions(j,1) && Stack_xy_aligned_positions(j,1)+sizeStack(j,1)>=Stack_xy_aligned_positions(i,1) && Stack_xy_aligned_positions(i,2)+sizeStack(i,2)>=Stack_xy_aligned_positions(j,2) && Stack_xy_aligned_positions(j,2)+sizeStack(j,2)>=Stack_xy_aligned_positions(i,2)
            Ctemp=zeros(1,sizeStack(i,3)+sizeStack(j,3)-1);
            Ntemp=zeros(size(Ctemp));
            nn=0;
            for k=-sizeStack(i,3)+1:sizeStack(j,3)-1
                nn=nn+1;
                [Ctemp(nn),Ntemp(nn)] = Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),[Stack_xy_aligned_positions(i,:),k+1],[Stack_xy_aligned_positions(j,:),1]);
            end
            [Cmax(i,j),Cind]=max(Ctemp);
            Nmax(i,j)=Ntemp(Cind);
            delZmax(i,j)=-sizeStack(i,3)+Cind;
        end
    end
end

Cmax=Cmax+Cmax';
Nmax=Nmax+Nmax';
delZmax=delZmax-delZmax';
W=Cmax.*Nmax.^0.5; % weight

Stack_z_aligned_positions = round((W-diag(sum(W)-10^-6))^(-1)*diag(W*delZmax));
Stack_z_aligned_positions=Stack_z_aligned_positions-min(Stack_z_aligned_positions)+1;
% temp=Stack_z_aligned_positions*ones(1,length(Stack_z_aligned_positions))-ones(length(Stack_z_aligned_positions),1)*Stack_z_aligned_positions';
% plot(temp(delZmax(:)>0),delZmax(delZmax(:)>0),'*')
% xlabel('optimal \Delta z'), ylabel('max \Delta z')


% align stacks in xyz
Stack_xyz_aligned_positions=[Stack_xy_aligned_positions,Stack_z_aligned_positions];
for k=1:0
    C=zeros(N_stacks,N_stacks);
    dCdxi=zeros(N_stacks,N_stacks);
    dCdyi=zeros(N_stacks,N_stacks);
    dCdzi=zeros(N_stacks,N_stacks);
    for i=1:N_stacks
        for j=(i+1):N_stacks
            C(i,j) = Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),Stack_xyz_aligned_positions(i,:),Stack_xyz_aligned_positions(j,:));
            dCdxi(i,j) = (Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),Stack_xyz_aligned_positions(i,:)+[1,0,0],Stack_xyz_aligned_positions(j,:))-Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),Stack_xyz_aligned_positions(i,:)+[-1,0,0],Stack_xyz_aligned_positions(j,:)))/2;
            dCdyi(i,j) = (Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),Stack_xyz_aligned_positions(i,:)+[0,1,0],Stack_xyz_aligned_positions(j,:))-Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),Stack_xyz_aligned_positions(i,:)+[0,-1,0],Stack_xyz_aligned_positions(j,:)))/2;
            dCdzi(i,j) = (Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),Stack_xyz_aligned_positions(i,:)+[0,0,1],Stack_xyz_aligned_positions(j,:))-Corr_Sparse(Stack{i},Stack{j},sizeStack(i,:),sizeStack(j,:),Stack_xyz_aligned_positions(i,:)+[0,0,-1],Stack_xyz_aligned_positions(j,:)))/2;
        end
    end
    Stack_xyz_aligned_positions=Stack_xyz_aligned_positions+sign([sum((dCdxi-dCdxi'),2),sum((dCdyi-dCdyi'),2),sum((dCdzi-dCdzi'),2)]);
    [k,mean(C(:))]
end


% combine stacks and desplay max projection
[Orig,sizeOrig]=Combine_Stacks(Stack,sizeStack,Stack_xyz_aligned_positions);
figure
PlotMaxProjection_Sparse([],Orig,sizeOrig)
hold on
title('xyz Aligned Mosaic')
for i=1:N_stacks
    text(Stack_xyz_aligned_positions(i,2)-min(Stack_xyz_aligned_positions(:,2))+1,Stack_xyz_aligned_positions(i,1)-min(Stack_xyz_aligned_positions(:,1))+1,num2str(i),'Color','r')
end

% Aligned Start Points
SP_stack_number=zeros(size(SP_xy,1),1);
SP_z=zeros(size(SP_xy,1),1);
for i=1:size(SP_xy,1)
    SP_stack_number(i)=find(SP_xy(i,1)>=Stack_xy_positions(:,1) & SP_xy(i,1)<Stack_xy_positions(:,1)+sizeStack(:,1) & SP_xy(i,2)>=Stack_xy_positions(:,2) & SP_xy(i,2)<Stack_xy_positions(:,2)+sizeStack(:,2),1,'first');
    [~,SP_z(i)]=max(Stack{SP_stack_number(i)}(sub2ind(sizeStack(SP_stack_number(i),1:2),SP_xy(i,1)-Stack_xy_positions(SP_stack_number(i),1)+1,SP_xy(i,2)-Stack_xy_positions(SP_stack_number(i),2)+1)+prod(sizeStack(SP_stack_number(i),1:2)).*[0:sizeStack(SP_stack_number(i),3)-1]));
end
SP_xyz_aligned=[SP_xy+Stack_xyz_aligned_positions(SP_stack_number,1:2)-Stack_xy_positions(SP_stack_number,:),SP_z+Stack_xyz_aligned_positions(SP_stack_number,3)-1];
plot(SP_xyz_aligned(:,2),SP_xyz_aligned(:,1),'r*')
%save([pth,'xyz_Aligned_Mosaic_3Median_0_005.mat'],'Orig','sizeOrig','SP_xyz_aligned','-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function combines multiple stacks in a single sparse stack
function [Im,sizeIm]=Combine_Stacks(Stack,sizeStack,Stack_xyz_positions)

N_stacks=length(Stack);
sizeIm=[max(Stack_xyz_positions+sizeStack)-min(Stack_xyz_positions)];
N_voxels=0;
for i=1:N_stacks
    N_voxels=N_voxels+nnz(Stack{i});
end
IND=zeros(N_voxels,1);
V=zeros(N_voxels,1);
count=0;

for i=1:N_stacks
    [ind,~,v]=find(Stack{i});
    [x,y,z]=ind2sub_AS(sizeStack,ind);
    
    x=x+Stack_xyz_positions(i,1)-min(Stack_xyz_positions(:,1));
    y=y+Stack_xyz_positions(i,2)-min(Stack_xyz_positions(:,2));
    z=z+Stack_xyz_positions(i,3)-min(Stack_xyz_positions(:,3));
    IND(count+1:count+length(v))=sub2ind(sizeIm,x,y,z);
    V(count+1:count+length(v))=v;
    count=count+length(v);
end
V(IND==0)=[];
IND(IND==0)=[];
temp=sortrows([IND,V],[1 2]);
[IND,m,~]=unique(temp(:,1),'last');
V=temp(m,2);
Im=sparse(IND,ones(length(IND),1),V);


% this function creates a sparse max projection of a stack
function MaxProjIm = MaxProjection(Im,sizeIm)

MaxProjIm=sparse(sizeIm(1),sizeIm(2));
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
    MaxProjIm(ind_ij)=temp_max;
end


% This function calculates normalized correlation, C_xy, between 2 sparse 
% projected stacks. N_xy is the number of overlaping pixels.
function [C_xy,N_xy] = Corr_Sparse_xy(StackProj1,StackProj2,size1,size2,R1,R2,thr)

C_xy=0; N_xy=0;

if R1(1)+size1(1)>=R2(1) && R2(1)+size2(1)>=R1(1) && R1(2)+size1(2)>=R2(2) && R2(2)+size2(2)>=R1(2)
    Overlap_min=max([R1;R2]);
    Overlap_max=min([R1+size1;R2+size2])-1;
    Overlap_size=Overlap_max-Overlap_min+1;
    
    StackProj1=StackProj1(Overlap_min(1)-R1(1)+1:Overlap_max(1)-R1(1)+1,Overlap_min(2)-R1(2)+1:Overlap_max(2)-R1(2)+1);
    StackProj2=StackProj2(Overlap_min(1)-R2(1)+1:Overlap_max(1)-R2(1)+1,Overlap_min(2)-R2(2)+1:Overlap_max(2)-R2(2)+1);
    StackProj1(StackProj1<thr)=0;
    StackProj2(StackProj2<thr)=0;
    
    if mean(StackProj1(:))*mean(StackProj2(:))~=0
        %C_xy=full((mean(mean(StackProj1.*StackProj2))));
        C_xy=full((mean(mean(StackProj1.*StackProj2))-mean(StackProj1(:))*mean(StackProj2(:)))./(var(StackProj1(:))*var(StackProj2(:)))^0.5);
        N_xy=nnz(StackProj1.*StackProj2);
    end
end


% This function calculates normalized correlation, C, between 2 sparse
% stacks. N is the number of overlaping voxels.
function [C,N] = Corr_Sparse(Stack1,Stack2,size1,size2,R1,R2)

C=0; N=0;

if R1(1)+size1(1)>=R2(1) && R2(1)+size2(1)>=R1(1) && R1(2)+size1(2)>=R2(2) && R2(2)+size2(2)>=R1(2) && R1(3)+size1(3)>=R2(3) && R2(3)+size2(3)>=R1(3)
    [ind1,~,v1]=find(Stack1);
    [x1,y1,z1]=ind2sub_AS(size1,ind1);
    [ind2,~,v2]=find(Stack2);
    [x2,y2,z2]=ind2sub_AS(size2,ind2);
    clear ind1 ind2
    
    Overlap_min=max([R1;R2]);
    Overlap_max=min([R1+size1;R2+size2])-1;
    Overlap_size=Overlap_max-Overlap_min+1;
    
    x1=x1-Overlap_min(1)+R1(1); y1=y1-Overlap_min(2)+R1(2); z1=z1-Overlap_min(3)+R1(3);
    x2=x2-Overlap_min(1)+R2(1); y2=y2-Overlap_min(2)+R2(2); z2=z2-Overlap_min(3)+R2(3);
    
    ind1_overlap=(x1>=1 & x1<=Overlap_size(1) & y1>=1 & y1<=Overlap_size(2) & z1>=1 & z1<=Overlap_size(3));
    ind2_overlap=(x2>=1 & x2<=Overlap_size(1) & y2>=1 & y2<=Overlap_size(2) & z2>=1 & z2<=Overlap_size(3));  
    
    x1=x1(ind1_overlap); y1=y1(ind1_overlap); z1=z1(ind1_overlap); v1=v1(ind1_overlap);
    x2=x2(ind2_overlap); y2=y2(ind2_overlap); z2=z2(ind2_overlap); v2=v2(ind2_overlap);
    clear ind1_overlap ind2_overlap
    
    temp1=sparse(sub2ind_AS(Overlap_size,x1,y1,z1),ones(length(v1),1),v1,prod(Overlap_size),1);
    temp2=sparse(sub2ind_AS(Overlap_size,x2,y2,z2),ones(length(v2),1),v2,prod(Overlap_size),1);
    if mean(temp1(:))*mean(temp2(:))~=0
        C=full((mean(temp1(:).*temp2(:))-mean(temp1(:))*mean(temp2(:)))./(var(temp1(:))*var(temp2(:)))^0.5);
        N=nnz(temp1.*temp2);
    end
end


