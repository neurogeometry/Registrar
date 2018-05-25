cluster = parpool('local',4);
inMatlab = mapreducer(0);
inPool = mapreducer(cluster);

ds = datastore('airlinesmall.csv','TreatAsMissing','NA',...
     'SelectedVariableNames','ArrDelay','ReadSize',1000);
preview(ds)

outputFolder = 'E:\temp';

meanDelay = mapreduce(ds,@meanArrivalDelayMapper,@meanArrivalDelayReducer,inMatlab);
readall(meanDelay)

meanDelay = mapreduce(ds,@meanArrivalDelayMapper,@meanArrivalDelayReducer,inPool,'OutputFolder',outputFolder);
readall(meanDelay)

delete(gcp)