% This function performs non-rigid stack and trace alignment
% It is based on a triangular elastic net

function [r1a,r2a]=TraceAlign3D(AMlbl1,r1,I1,AMlbl2,r2,I2)

r1a=[]; r2a=[];

k_mesh=1;
pad=50; % in voxels
MeshSize=5; % in voxels
betta=0.01;

%%%% create a triangular mesh
b1=[1,0,0]'; b2=[1/2,3^0.5/2,0]'; b3=[1/2,3^0.5/6,(2/3)^0.5]';
B=[b1,b2,b3].*MeshSize;
bc=[1,3^0.5/2,(2/3)^0.5].*MeshSize;
Range_min=min([min(r1);min(r2)])-pad;
Range_max=max([max(r1);max(r2)])+pad;
Range_min(1)=Range_min(1)-(Range_max(1)-Range_min(1))/3^0.5;
N=ceil((Range_max-Range_min)./bc);

ind_mesh=(1:prod(N))';
[x,y,z]=ind2sub(N,ind_mesh);
r_mesh=([x,y,z]-1)*B'+ones(size(x))*Range_min;
Neighb_xyz=[1,0,0;0,1,0;0,0,1;1,-1,0;1,0,-1;0,1,-1;];
ii=zeros(6*length(ind_mesh),1);
jj=zeros(6*length(ind_mesh),1);
count=1;
for i=1:length(ind_mesh)
    xtemp=x(i)+Neighb_xyz(:,1);
    ytemp=y(i)+Neighb_xyz(:,2);
    ztemp=z(i)+Neighb_xyz(:,3);
    indtemp=xtemp+(ytemp-1).*N(1)+(ztemp-1).*(N(1)*N(2));
    indtemp=indtemp(xtemp>=1 & xtemp<=N(1) & ytemp>=1 & ytemp<=N(2) & ztemp>=1 & ztemp<=N(3));
    ii(count:count+length(indtemp)-1)=ones(length(indtemp),1).*i;
    jj(count:count+length(indtemp)-1)=indtemp;
    count=count+length(indtemp);
end
ii=ii(1:count-1);
jj=jj(1:count-1);
AM_mesh=sparse(ii,jj,ones(1,length(i)),length(ind_mesh),length(ind_mesh));
AM_mesh=max(AM_mesh,AM_mesh');
%Nedges=sum(AM_mesh(:))/2;
L_mesh=AM_mesh-diag(sum(AM_mesh));

%%%% Define Boundary, BC
BC=sum(AM_mesh)<12;

%%%% find barycentric coordinates of all trace vertices
tr1=[(r1(:,1)-Range_min(1))./bc(:,1),(r1(:,2)-Range_min(2))./bc(:,2),(r1(:,3)-Range_min(3))./bc(:,3)];
t1=r1*(B*(B'*B)^-1); %[t1,t2,t3]


%%%% Gradient Descent
r_mesh(~BC,:)=r_mesh(~BC,:)+MeshSize.*rand(size(r_mesh(~BC,:)));
U_mesh0=-(k_mesh/2).*(r_mesh(:,1)'*L_mesh*r_mesh(:,1)+r_mesh(:,2)'*L_mesh*r_mesh(:,2)+r_mesh(:,3)'*L_mesh*r_mesh(:,3)); 
U_mesh=U_mesh0;
delU_mesh=Inf;
Nsteps=0;
while abs(delU_mesh)>0.000000000001*U_mesh0 && Nsteps<=10000
    Nsteps=Nsteps+1;
    disp([Nsteps,U_mesh])
    delr_mesh=betta.*L_mesh*r_mesh;
    delr_mesh(BC,:)=0;
    r_mesh=r_mesh+delr_mesh;
    U_mesh_old=U_mesh;
    U_mesh=-(k_mesh/2).*(r_mesh(:,1)'*L_mesh*r_mesh(:,1)+r_mesh(:,2)'*L_mesh*r_mesh(:,2)+r_mesh(:,3)'*L_mesh*r_mesh(:,3)); 
    delU_mesh=U_mesh-U_mesh_old;
end


% figure
% PlotAM(AM_mesh,r_mesh), hold on
% plot3(r_mesh(:,2),r_mesh(:,1),r_mesh(:,3),'b.')
% %xlim([0 10]), ylim([0 10]), zlim([0 10])
% axis equal, box on
% PlotAM(AMlbl1,r1)
% PlotAM(AMlbl2,r2)



