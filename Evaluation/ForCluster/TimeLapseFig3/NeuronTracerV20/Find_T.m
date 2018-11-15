% this function finds parameter T used in the confidence measure
% Costs is a cell array of merger costs within each cluster.
% Ec is the target cluster error-rate

function T=Find_T(Costs,Ec)

Tg=(0.01:0.01:5)';

F=zeros(size(Tg));
for i=1:length(Costs)
    F=F+exp(-min(Costs{i})./Tg)./sum(exp(-(1./Tg)*Costs{i}),2);
end

F=1-F/length(Costs);
[~,ind]=min(abs(F-Ec));
T=Tg(ind);
