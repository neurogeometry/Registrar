function [H,d] = Transformation3D(TargetLocations,SourceLocations,RandomSamples)

RandTargetLocations = TargetLocations(:,RandomSamples);
RandSourceLocations = SourceLocations(:,RandomSamples);
[M,N] = size(RandTargetLocations);
TargetExpand = zeros(N*M,M*M);
K=0;
q = 1;
for i = 1:N
    for j = 1:M
          TargetExpand(q,K*M+1:(K*M+M)) = RandTargetLocations(:,i)';
          K=K+1;
          q = q+1;
    end
    K=0;
end
H = reshape(TargetExpand\RandSourceLocations(:), M, M)';

 TransformedPoints = H*TargetLocations;
 d = sum((SourceLocations-TransformedPoints).^2,1);
end