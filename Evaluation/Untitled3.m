load('E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\result_NR_fiji.mat')
mean_dist = zeros(12,3);
figure
for sourceID = 1:12
    mean_dist(sourceID,1) = mean(mean((result{sourceID}.Bouton.r2-result{sourceID}.Bouton.r1).^2,1).^0.5);
    mean_dist(sourceID,2) = mean(mean((result{sourceID}.Bouton.r2_fiji_Affine-result{sourceID}.Bouton.r1_fiji_Affine).^2,1).^0.5);
    mean_dist(sourceID,3) = mean(mean((result{sourceID}.Bouton.r2-result{sourceID}.Bouton.r1_NR).^2,1).^0.5);
end
plot(mean_dist(:,1),'o')
hold on
plot(mean_dist(:,2),'x')
hold on
plot(mean_dist(:,3),'s')