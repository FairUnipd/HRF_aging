clear all
clc
close all

base_path=pwd;
files=dir(fullfile(base_path,'FINAL_RESULT_folder'));
addpath(genpath(files(1).folder))
load node_lab_conc.mat

labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];

%load SUM_COV generated for Figure6
load multivar_sel.mat


%% select top 5 features
[feat_top,idx_top]=sort(SUM_cov,'descend');
selec_final_feat=sort(idx_top(1:5));

%% retrain minimal model

%load initial X, Y and Y_LR with all subjects
load DATA_LOGREG.mat

for ss=1:1000
    
    %set random seed
    rng(ss,'Threefry');
    
    partition = cvpartition(Y,'Holdout',0.2,'Stratify',true);
    Train_idx = training(partition);
    Test_idx = test(partition);
    
    % split the input
    XTest = X(Test_idx,selec_final_feat);
    XTrain=X(Train_idx,selec_final_feat);
    % split the output
    YTest = Y(Test_idx,:);
    YTrain=Y(Train_idx,:);
    Y_LRTrain=Y_LR(Train_idx,:);

    XTrain = (XTrain-mean(XTrain))./std(XTrain);
    XTest  = (XTest-mean(XTrain))./std(XTrain);
    
    %glm model
    varNames_S = matlab.lang.makeValidName([CovNames(selec_final_feat),{'Output'}]);
    minimal_model=fitglm(XTrain, YTrain, 'linear', 'Distribution', 'binomial',...
            'link', 'logit','dispersionFlag',true,'VarNames',varNames_S');
    
    b_S=minimal_model.Coefficients.Estimate; 
    p_S=glmval(b_S,XTest,'logit');

    %model goodness metrics
    zscore_test(ss,:)=minimal_model.Coefficients.pValue;
    beta_test(ss,:)=minimal_model.Coefficients.Estimate;
    AUC_test(ss)=c_index(p_S,YTest);
    
    %measure performance metrics
     [~, ~, ~, auc_roc,auc_prc, ~, balanced_accuracy_opt_ba,...
    ~, ~,~,accuracy_opt_ba] = performance_wBA(p_S,YTest);
    
 
    AUC_ROC(ss)=auc_roc;
    AUC_PRC(ss)=auc_prc;
    BA_ROC(ss)=balanced_accuracy_opt_ba;
    ACC_ROC(ss)=accuracy_opt_ba;
  

    
end
%clear vars
clear x_roc y_roc y_prc auc_roc auc_prc cutoff_opt_ba balanced_accuracy_opt_ba...
    sensitivity_opt_ba specificity_opt_ba precision_opt_ba accuracy_opt_ba

%retain and combine only good models
beta_test(zscore_test>=0.05)=nan;
beta_test(AUC_test<=0.5,:)=nan;
sum(isnan(beta_test))

%% Test minimal model on Dataset 1

XTest = X(:,selec_final_feat);
YTest_all = Y;

%variable normalization
sd_lem=std(XTest,'omitnan');
mean_lem=mean(XTest,'omitnan');
mean_balanced_lem=mean([mean(XTest(1:sum(Y==1),:));mean(XTest(sum(Y==1)+1:end,:))]);
XTest  = (XTest-mean_balanced_lem)./std(XTest); %zscore

%get median beta coefficients and test on Dataset 1
b_S=median(beta_test,'omitnan')'; 
p_S_all=glmval(b_S,XTest,'logit');
[x_roc, y_roc, y_prc, auc_roc,auc_prc, cutoff_opt_ba, balanced_accuracy_opt_ba,...
    sensitivity_opt_ba, specificity_opt_ba,precision_opt_ba,accuracy_opt_ba] = performance_wBA(p_S_all,YTest_all);
    

%% Boxplot of beta coefficients

load COVNAMES_struct.mat
load group_scatter.mat

for jj=1:length(selec_final_feat)

    ROIs_sel(jj)=CovNames_struct(selec_final_feat(jj)).ROI;
    HRF_sel=CovNames_struct(selec_final_feat(jj)).HRF_param;
    HRF_params_selec(jj)=find(matches(labels_hrf,HRF_sel));

end

cortical_ROIS=find(matches(node_labels,ROIs_sel));
cortical_ROIS=cortical_ROIS(1:end-1);

%RSN color coding
group_col = [0.9608    0.9294    0.9294;0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
    0.3922    0.8314    0.0745;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
    0.8392    0.5922    0.0157;0.850980392156863   0.325490196078431   0.098039215686275];
marker_list={'o','+','x','s','^','diamond','v','pentagram','>'};
color_list=gray(9);
reorder=[3 5 2 1 4];

beta_plot=beta_test(:,2:end);
beta_plot=beta_plot(:,reorder);
HRF_params_selec_reorder=HRF_params_selec(reorder);
col_reordered=[group_col(4,:);[0 1 0.5];group_col(4,:);group_col(2,:);group_col(7,:)]


figure
hold on
b=boxplot(beta_plot,'Colors',col_reordered,'Widths',0.5)
set(b,{'linew'},{1.5})
hold off
a= get(get(gca,'children'),'children');
for mm=1:length(ROIs_sel)
    

    a(abs(6-mm)).Marker=marker_list(HRF_params_selec_reorder(mm));
    a(abs(6-mm)).MarkerEdgeColor=[0.380392156862745   0.352941176470588   0.352941176470588];
    a(abs(6-mm)).MarkerFaceColor=color_list(HRF_params_selec_reorder(mm),:);
    a(abs(6-mm)).LineWidth=1;
    
end

ylim([-10 10])
xticks(1:5)
xticklabels(ROIs_sel(reorder))

%% project onto surfaces

%cortical surfaces
addpath(genpath('ENIGMA'))

for rr=1:length(cortical_ROIS)
    HRF_selection_sum_cortical=zeros(62,1)';
    HRF_selection_sum_cortical(cortical_ROIS(rr))=1;
    group_num_filtered=group_num(1:62).*(HRF_selection_sum_cortical>0);
    corr_map_conte69 = parcel_to_surface(group_num_filtered,'schaefer_62_conte69');
    col_map=group_col([1:nonzeros(group_num_filtered)+1],:);
    % Project the results on the surface brain
    g= figure,plot_cortical_col_map(corr_map_conte69,col_map, 'surface_name', 'conte69')
end

%subcortical surfaces
HRF_selection_sum_subcortical=zeros(10,1)';
HRF_selection_sum_subcortical(6)=1;
reorder=[2 5 4 3 1 7 10 9 8 6];
HRF_selection_sum_subcortical=HRF_selection_sum_subcortical(reorder);
group_num=[1 2 3 4 5 1 2 3 4 5];
group_num_filtered=group_num.*(HRF_selection_sum_subcortical>0);
g_sub=figure,plot_subcortical_AAL3(group_num_filtered,'ventricles','False')


%% PDFs plots

col=lines(4);

figure
h1=histogram(AUC_ROC,'BinWidth',0.04,'FaceAlpha',0.2,'FaceColor',col(1,:),'Normalization','pdf','BinLimits',[0 1]);
hold on
h1_curve=interp1([h1.BinEdges(1:end-1)+h1.BinWidth/2 1+10e-3],[h1.Values 0],h1.BinLimits(1):0.01:h1.BinLimits(2)+10e-3,'spline','extrap');
h1_bounds=h1.BinLimits(1):0.01:h1.BinLimits(2)+10e-3;
h1_curve(h1_curve<0)=0;

plot(h1.BinLimits(1):0.01:h1.BinLimits(2)+10e-3,h1_curve);

figure
h2=histogram(AUC_PRC,'BinWidth',0.04,'FaceAlpha',0.2,'FaceColor',col(2,:),'Normalization','pdf','BinLimits',[0 1]);
hold on
h2_curve=interp1([h2.BinEdges(1:end-1)+h2.BinWidth/2 1+10e-3],[h2.Values 0],h2.BinLimits(1):0.01:h2.BinLimits(2)+10e-3,'spline','extrap');
h2_bounds=h2.BinLimits(1):0.01:h2.BinLimits(2)+10e-3;
h2_curve(h2_curve<0)=0;
plot(h2.BinLimits(1):0.01:h2.BinLimits(2)+10e-3,h2_curve);

figure
hold on
h3=histogram(BA_ROC,'BinWidth',0.04,'FaceAlpha',0.2,'FaceColor',col(3,:),'Normalization','pdf','BinLimits',[0 1]);
hold on
h3_curve=interp1([h3.BinEdges(1:end-1)+h3.BinWidth/2 1+10e-3],[h3.Values 0],h3.BinLimits(1):0.01:h3.BinLimits(2)+10e-3,'spline','extrap');
h3_bounds=h3.BinLimits(1):0.01:h3.BinLimits(2)+10e-3
h3_curve(h3_curve<0)=0;
plot(h3.BinLimits(1):0.01:h3.BinLimits(2)+10e-3,h3_curve);

figure
hold on
h4=histogram(ACC_ROC,'BinWidth',0.04,'FaceAlpha',0.2,'FaceColor',col(4,:),'Normalization','pdf','BinLimits',[0 1]);
hold on
h4_curve=interp1([h4.BinEdges(1:end-1)+h4.BinWidth/2 1+10e-3],[h4.Values 0],h4.BinLimits(1):0.01:h4.BinLimits(2)+10e-3,'spline','extrap');
h4_curve(h4_curve<0)=0;
h4_bounds=h4.BinLimits(1):0.01:h4.BinLimits(2)+10e-3;
plot(h4.BinLimits(1):0.01:h4.BinLimits(2)+10e-3,h4_curve);


%plot shaded curves
figure
plot_shaded(h1_bounds,h1_curve,'Color',col(1,:),'LineWidth',2)
hold on
plot_shaded(h2_bounds,h2_curve,'Color',col(2,:),'LineWidth',2)
hold on
plot_shaded(h3_bounds,h3_curve,'Color',col(3,:),'LineWidth',2)
hold on
plot_shaded(h4_bounds,h4_curve,'Color',col(4,:),'LineWidth',2)
xlim([0.2 1])
ylabel('PDF')

%% performance metrics on Dataset 1 (for Figure 8)

%C-index
AUC_test=c_index(p_S_all,YTest_all)

%confusion chart
yhat_S_all = (p_S_all>=cutoff_opt_ba); %cutoff_opt_roc_S/cutoff_opt_prc_S
conf_mat=confusionmat(logical(YTest_all),yhat_S_all);
figure
c_S = confusionchart(conf_mat,categorical({'Old','Young'}));
balanced_acc=((conf_mat(1,1)/sum(conf_mat(1,:)))+(conf_mat(2,2)/sum(conf_mat(2,:))))/2
c_S.ColumnSummary = 'column-normalized';
c_S.RowSummary = 'row-normalized';
c_S.Title=['CONFUSION MATRIX- Balanced accuracy=' num2str(balanced_acc)];
c_S.DiagonalColor=[0.074509803921569   0.623529411764706   1];
c_S.OffDiagonalColor=[1.0000    1 0.066666666666667];

%ROC curve
figure(81)
p1=plot(x_roc,y_roc,'-','linewidth',2,'color',[0.074509803921569   0.623529411764706   1])
hold on
plot([0:1],[0:1],'k--')
hold on
p2=plot(1-specificity_opt_ba,sensitivity_opt_ba,'bo','markersize',10,'markerfacecolor','b','linewidth',2)
set(gca,'fontsize',12)
sgtitle(['ROC curve-AUC=',num2str(AUC_test)])
xlabel('1-specificity')
ylabel('sensitivity or recall')
hold on
legend([p1,p2],{'ROC Train set','Random classifier'})

%PR curve
P=sum(YTest_all==1);
N=sum(YTest_all==0);
figure(82)
p3=plot(y_roc,y_prc,'-','linewidth',2,'color',[0.074509803921569   0.623529411764706   1])
hold on
plot(sensitivity_opt_ba,precision_opt_ba,'bo','markersize',10,'markerfacecolor','b','linewidth',2)
set(gca,'fontsize',12)
sgtitle('PR curve')
xlabel('sensitivity or recall')
ylabel('precision')
hold on
p4=yline(P/(P+N),'--','Color',[0.074509803921569   0.623529411764706   1])
legend([p3 p4],{'PRC Train set','Random classifier Train set'})


