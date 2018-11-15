% This function calculates the average distance between two setts of traces
% The two setts should contain equal numbers of traces

function [E_um, E_voxel] = TraceDistance(AMg, rg, AMa, ra,pixelSize,ShowPlot)
%Distances = [];
UnassignedWorkerCost=5;
UnassignedJobCost=5;
    
%relabel AMa based on distance
Lg=unique(AMg(AMg>0));
E_voxel=nan(1,length(Lg));
E_um=nan(1,length(Lg));
for i=1:length(Lg)
    [e1,e2]=find(AMg==Lg(i));
    ee=unique([e1,e2]);
    AMg_temp=AMg(ee,ee);
    rg_temp=rg(ee,:);
      Dag_temp=(bsxfun(@minus,ra(:,1),rg_temp(:,1)').^2+bsxfun(@minus,ra(:,2),rg_temp(:,2)').^2+bsxfun(@minus,ra(:,3),rg_temp(:,3)').^2).^0.5;
%     Dag_temp=(((bsxfun(@minus,ra(:,1),rg_temp(:,1)'))*pixelSize(1)).^2+((bsxfun(@minus,ra(:,2),rg_temp(:,2)'))*pixelSize(2)).^2+((bsxfun(@minus,ra(:,3),rg_temp(:,3)'))*pixelSize(3)).^2).^0.5;
    
    
    [~,ind]=min(Dag_temp(:));
    [f,~]=ind2sub(size(Dag_temp),ind);
    L_temp=max(AMa(f,:));
    
    [e1,e2]=find(AMa==L_temp);
    ee=unique([e1,e2]);
    AMa_temp=AMa(ee,ee);
    ra_temp=ra(ee,:);
    Dag_temp=Dag_temp(ee,:);
    
    [AM,~]=Hungarian_fast(Dag_temp,UnassignedWorkerCost,UnassignedJobCost);
    %Distances = [Distances;Dag_temp(AM)];
    [ee1,ee2]=find(AM);
    E_um(i)=mean(sum(((ra_temp(ee1,:)-rg_temp(ee2,:)).*(ones(length(ee1),1)*pixelSize)).^2,2).^0.5);
    E_voxel(i)=mean(sum((ra_temp(ee1,:)-rg_temp(ee2,:)).^2,2).^0.5);
    %E_voxel(i)=mean(Dag_temp(AM));
    
    if isnan(mean(Dag_temp(AM)))
        while isnan(mean(Dag_temp(AM)))
            UnassignedWorkerCost=UnassignedWorkerCost+1;
            UnassignedJobCost=UnassignedJobCost+1;
            [AM,~]=Hungarian_fast(Dag_temp,UnassignedWorkerCost,UnassignedJobCost);
            %Distances = [Distances;Dag_temp(AM)];
%             E_voxel(i)=mean(Dag_temp(AM));
            [ee1,ee2]=find(AM);
            E_um(i)=mean(sum(((ra_temp(ee1,:)-rg_temp(ee2,:)).*(ones(length(ee1),1)*pixelSize)).^2,2).^0.5);
            E_voxel(i)=mean(sum((ra_temp(ee1,:)-rg_temp(ee2,:)).^2,2).^0.5);
        end
    end
    
    if ShowPlot
        figure
        PlotAM(AMg_temp,rg_temp,'r')
        axis equal
        PlotAM(AMa_temp,ra_temp,'g')
        axis equal
        
        [e1,e2]=find(AM);
        X=[ra_temp(e1,1),rg_temp(e2,1)]';
        Y=[ra_temp(e1,2),rg_temp(e2,2)]';
        Z=[ra_temp(e1,3),rg_temp(e2,3)]';
    line(Y,X,Z,'Color','k','LineWidth',1);
    end
end



