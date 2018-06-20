function FeatureExtractionVisualization(path,pints)

IM_Original=ImportStack(path);
figure,imshow(max(IM_Original,[],3),[0 max(IM_Original(:))])
hold on
plot(pints(:,2),pints(:,1),'r*');


end