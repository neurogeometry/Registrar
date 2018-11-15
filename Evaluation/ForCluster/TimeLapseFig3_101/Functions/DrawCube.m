function DrawCube(origin, boxsize)
 
%  vec = 'ymcrgbk';
%  clr = vec(randi(numel(vec),N));

%     x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*boxsize(1)+origin(1)+boxsize(1)/2;
%     y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*boxsize(2)+origin(2)+boxsize(2)/2;
%     z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*boxsize(3)+origin(3)+boxsize(3)/2;

x=[0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0]*boxsize(1)+origin(1);
y=[0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1]*boxsize(2)+origin(2);
z=[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1]*boxsize(3)+origin(3);
h=line([x(1:2:end);x(2:2:end)],[y(1:2:end);y(2:2:end)],[z(1:2:end);z(2:2:end)],'Color','b');
%     for j=1:6
%         h=line(x(:,j),y(:,j),z(:,j));
% %          h=patch(x(:,j),y(:,j),z(:,j),'w');
%           set(h,'edgecolor','k')
% %          set(h,'facecolor',clr(i))
%     end
end