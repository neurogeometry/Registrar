function DrawStackMap(StackPositions_pixels, StackSizes_pixels)
 
x=StackSizes_pixels(:,1)*[0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0]+StackPositions_pixels(:,1)*ones(1,24);
y=StackSizes_pixels(:,2)*[0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1]+StackPositions_pixels(:,2)*ones(1,24);
z=StackSizes_pixels(:,3)*[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1]+StackPositions_pixels(:,3)*ones(1,24);

xstart=x(:,1:2:end-1);
xend=x(:,2:2:end);
ystart=y(:,1:2:end-1);
yend=y(:,2:2:end);
zstart=z(:,1:2:end-1);
zend=z(:,2:2:end);

line([xstart(:)';xend(:)'],[ystart(:)';yend(:)'],[zstart(:)';zend(:)'],'Color','b');
axis equal
drawnow
hold on
%     [~,StackName]=fileparts(StackList{j,1});
%     text(StackPositions_pixels(j,1)+boxsize(1)/2,StackPositions_pixels(j,2)+boxsize(2)/2,StackPositions_pixels(j,3)+boxsize(3)/2,[num2str(j),'-',StackName],'FontSize',6),hold on



end