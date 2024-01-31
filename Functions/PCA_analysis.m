function [tot_perc,pc1_young,pc2_young,pc3_young,pc1_old,pc2_old,pc3_old] = PCA_analysis(DATA,TR,HRF_length,hrf_young,hrf_old)
%performs and save results of PCA

[coeff,score,~,~,explained,~] = pca(DATA); 

coeff_pc1=coeff(:,1);
coeff_pc2=coeff(:,2);
coeff_pc3=coeff(:,3);

%display first 3 PCS
figure
TR=1.4;
F=0:TR:HRF_length*TR-TR;
plot(F,coeff_pc1,'k','LineWidth',2)
hold on
plot(F,coeff_pc2,'k--','LineWidth',2)
hold on
plot(F,coeff_pc3,'k-o','LineWidth',2)
xlabel('Time (s)')
legend('PC1','PC2','PC3')
hold off

%percentage of explained variance
perc_pc1=explained(1)
perc_pc2=explained(2)
perc_pc3=explained(3)
tot_perc=perc_pc1+perc_pc2+perc_pc3;

%PCA scores
pc1_young=score(1:size(hrf_young,2),1);
pc2_young=score(1:size(hrf_young,2),2);
pc3_young=score(1:size(hrf_young,2),3);

pc1_old=score(size(hrf_young,2)+1:size(hrf_young,2)+size(hrf_old,2),1);
pc2_old=score(size(hrf_young,2)+1:size(hrf_young,2)+size(hrf_old,2),2);
pc3_old=score(size(hrf_young,2)+1:size(hrf_young,2)+size(hrf_old,2),3);

end