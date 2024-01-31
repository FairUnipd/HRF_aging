clear all
clc

%figures to keep from previous code
figs2keep = [81,82];
all_figs = findall(0, 'type', 'figure');
all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));

%load info minimal model
load info_minimal_model.mat
load node_lab_conc.mat



labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];

%covariate labels
for ii=0:length(node_labels)-1
        for j=1:length(labels_hrf)
             CovNames(floor(ii*9+j))=strcat(node_labels(ii+1),'_',labels_hrf(j));
        end
end  
node_bounds=1:74;

%load HCP data (Dataset 2)
files=dir('RESULTS_POST_HCA_HCD_interp');
files=files(3:end);

%% create HRF partition

age_young=40;
age_old=60;
[params_y,params_o,HRFmatrix_y,HRFmatrix_o,Y_LR_y,Y_LR_o] = HRF_partition(files,age_young,age_old);


yo=reshape(params_y,666,size(Y_LR_y,2));
old=reshape(params_o,666,size(Y_LR_o,2));
X = [yo,old]';     %columns of covariates
Y= [ones(size(Y_LR_y,2),1); zeros(size(Y_LR_o,2),1)]; %columns of outputs
Y_LR=[Y_LR_y,Y_LR_o]';

XTest = X(:,selec_final_feat);

%variable normalization
sd_lem_forHCP=sd_lem;
sd_hcp=std(XTest,'omitnan');
mean_balanced=mean([mean(XTest(1:length(Y_LR_y),:));mean(XTest(length(Y_LR_y)+1:end,:))]);
XTest=(XTest-mean_balanced)./sd_hcp;   
YTest_all = Y;

%evaluate performance
p_S_all=glmval(b_S,XTest,'logit');

[x_roc, y_roc, y_prc, auc_roc,auc_prc, cutoff_opt_ba, balanced_accuracy_opt_ba,...
    sensitivity_opt_ba, specificity_opt_ba,precision_opt_ba,accuracy_opt_ba] = performance_wBA(p_S_all,YTest_all);
    


%% performance metrics on Dataset 2 


%C-index
AUC_test=c_index(p_S_all,YTest_all)

%confusion chart
yhat_S_all = (p_S_all>=cutoff_opt_ba); 
conf_mat=confusionmat(logical(YTest_all),yhat_S_all);
figure
c_S = confusionchart(conf_mat,categorical({'Old','Young'}));
balanced_acc=((conf_mat(1,1)/sum(conf_mat(1,:)))+(conf_mat(2,2)/sum(conf_mat(2,:))))/2
c_S.ColumnSummary = 'column-normalized';
c_S.RowSummary = 'row-normalized';
c_S.Title=['CONFUSION MATRIX- Balanced accuracy=' num2str(balanced_acc)];
c_S.DiagonalColor=[0.3922    0.8314    0.0745];
c_S.OffDiagonalColor=[1.0000    0.4118    0.1608];

%ROC curve
figure(111)
hold on
p1bis=plot(x_roc,y_roc,'-','linewidth',2,'color',[0.3922    0.8314    0.0745])
hold on
plot([0:1],[0:1],'k--')
hold on
plot(1-specificity_opt_ba,sensitivity_opt_ba,'bo','markersize',10,'markerfacecolor',[0.1294    0.6118    0.2902],'linewidth',2)
set(gca,'fontsize',12)
sgtitle(['ROC curve-AUC=',num2str(AUC_test)])
xlabel('1-specificity')
ylabel('sensitivity or recall')
xlim([-0.02 1.02])
ylim([-0.02 1.02])
hold on

%PR curve
P=sum(YTest_all==1);
N=sum(YTest_all==0);
figure(112)
hold on
p3bis=plot(y_roc,y_prc,'-','linewidth',2,'color',[0.3922    0.8314    0.0745])
hold on
plot(sensitivity_opt_ba,precision_opt_ba,'bo','markersize',10,'markerfacecolor',[0.1294    0.6118    0.2902],'linewidth',2)
set(gca,'fontsize',12)
sgtitle('PR curve')
xlabel('sensitivity or recall')
ylabel('precision')
xlim([-0.02 1.02])
ylim([0.4 1.02])
hold on
p4bis=yline(P/(P+N),'--','Color',[0.3922    0.8314    0.0745])

