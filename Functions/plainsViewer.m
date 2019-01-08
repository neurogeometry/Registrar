function plainsViewer(VisualizationStackHandle,IM)
%This function takes in a 3D image IM and 3-column position vector r. r can
%be empty. It allows the user to scroll through different planes, with
%points in the current plane highlighted.
%Examples:
%im3dscroll(rand(200,200,10),[(1:200)',(1:200)',(0.55:0.05:10.50)']);caxis([0 5]);
%im3dscroll(rand(200,200,10),[]);caxis([0 5]);
%Author:%Rohan, 1/3/2017, Matlab2016a.
% Seyed October 5, 2017
% hf=figure;
% hf.Visible='on';

% hf5 = findobj(NCT_Registration,'Tag', 'slider1');
% hf5 = VisualizationStackHandle.Children(1);
% hf5.Max = size(IM,3);
% hf5.Value = 1;
% hf5.Callback = @ScrollFcn;

% tb2 = findobj(NCT_Registration,'Tag', 'axes1');
tb2 = VisualizationStackHandle.Children(3);
hf = tb2;
% % tb2.Units='normalized';
h_im=imshow(IM(:,:,1),'Parent',tb2);hold on
hf.CLim = [0 max(IM(:))];
% hf.YLabel.String='X axis';hf.XLabel.String='Y axis';
hf.YLabel.Position(1) = 0;
hf.XLabel.Position(2) = -50;
hf.XLabel.Position(1) = 70;
pt2 = get(tb2,{'Position','tightinset','PlotBoxAspectRatio'});
pt2{3}(2);
%  hf.XLim=[0 900];hf.YLim=[0 1200];
tb2.Tag='axes1';


% tb3 = findobj(hf.Parent,'-depth',1,'Tag', 'axes5');
% hf3 = tb3;
% IM3 = squeeze(IM(:,:,:));
% h_im3=imshow(squeeze(IM(:,1,:)),'Parent',tb3);hold on
% tb3.Units='pixels';
% tb3.Tag='axes5';
% hf3.CLim = [0 max(IM3(:))];
% hf3.YLabel.String='X axis';hf3.XLabel.String='Z axis';
% hf3.XLim = [0.5 150];
% % hf3.YLim = [0.5 1100];
% hf3.YLabel.Position(1) = 0;
% hf3.XLabel.Position(2) = -50;
% 
% tb4 = findobj(hf.Parent,'-depth',1,'Tag', 'axes6');
% hf4 = tb4;
% IM4 = squeeze(IM(:,:,:));
% h_im4=imshow(imrotate(squeeze(IM(:,1,:)),90),'Parent',tb4);hold on
% tb4.Units='normalized';
% tb4.Tag='axes6';
% hf4.CLim = [0 max(IM4(:))];
% hf4.YLabel.String='Z axis';hf4.XLabel.String='Y axis';
% hf4.Units='pixels';%hf4.Position=[165, 50, 750, 100];
% hf4.XLim = [0 pt2{3}(1)*2];
% hf4.YLabel.Position(1) = 0;
% hf4.XLabel.Position(2) = -50;
% hf4.XLabel.Position(1) = 50;
% % hf4.YLim = hf.YLim; 

 
% set(hf4, 'Units', 'pixels', 'Position', [70, 100, 900, 50]);

% hf4.YLim=[0 150];

hf.UserData.IM=IM;
hf.UserData.currplane=1;
hf.UserData.h_im=h_im;

% hf3.UserData.IM=IM3;
% hf3.UserData.currplane=1;
% hf3.UserData.h_im=h_im3;
% 
% hf4.UserData.IM=IM4;
% hf4.UserData.currplane=1;
% hf4.UserData.h_im=h_im4;


tb2.Title.String=['Current plane: ',num2str(hf.UserData.currplane),' / ',num2str(size(hf.UserData.IM,3))];
ff=gcf;
tb2.Tag='axes1';
ff.WindowScrollWheelFcn=@ScrollFcn;
tb2.Tag='axes1';
%The figure has a property called WindowScrollWheelFcn that can be set to
%any function. A handle to such a function is passed, and the function is
%referred to as a callback function.

end

function ScrollFcn(src,ed)
%Callback function that is executed whenever the user scrolls. Callback
%functions that are executed when particular actions are performed (in this
%case a scroll operation) have 2 default inputs. Here the src is the figure
%handle, and ed is the event data. For WindowScrollWheelFcn, the default
%eventdata is a structure containing VerticalScrollCount along with other
%information.
% try
hf = gca;
% hf3 = findobj(hf.Parent,'-depth',1,'Tag', 'axes5');
% hf4 = findobj(hf.Parent,'-depth',1,'Tag', 'axes6');
% hf5 = findobj(hf.Parent,'-depth',1,'Tag', 'slider1');
% % hf.UserData.currplane = hf5.Value;

if ~isempty(hf.UserData)
    hi=hf.UserData.h_im;
%     hi3=hf3.UserData.h_im;
%     hi4=hf4.UserData.h_im;
    if ~isempty(hi)
        prevplane=hf.UserData.currplane;
        if strcmp(ed.EventName,'WindowScrollWheel')
            scrlnum = ed.VerticalScrollCount;
        else
            scrlnum = ed.Source.Value-hf.UserData.currplane;
        end
        if scrlnum>0
            hf.UserData.currplane=min(hf.UserData.currplane+scrlnum,size(hf.UserData.IM,3));
%             hf3.UserData.currplane=min(hf3.UserData.currplane+ed.VerticalScrollCount,size(hf3.UserData.IM,3));
%             hf4.UserData.currplane=min(hf4.UserData.currplane+ed.VerticalScrollCount,size(hf4.UserData.IM,3));
        else
            hf.UserData.currplane=max(hf.UserData.currplane+scrlnum,1);
%             hf3.UserData.currplane=max(hf3.UserData.currplane+ed.VerticalScrollCount,1);
%             hf4.UserData.currplane=max(hf4.UserData.currplane+ed.VerticalScrollCount,1);
            
        end
        hi.CData=hf.UserData.IM(:,:,round(hf.UserData.currplane));
%         hi3.CData=squeeze(hf.UserData.IM(:,hf.UserData.currplane,:));
%         hi4.CData=imrotate(squeeze(hf.UserData.IM(:,hf.UserData.currplane,:)),90);
%           hf5.Value = hf.UserData.currplane;
    end
    
    hf.Title.String=['Current plane: ',num2str(round(hf.UserData.currplane)),' / ',num2str(size(hf.UserData.IM,3))];
    
end
% catch
%     warndlg('Please click on the Image','!! Warning !!');
% end

end
