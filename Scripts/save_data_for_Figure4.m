clear all
clc
close all

load node_lab_conc.mat
labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];
node_bounds=1:74;

%load data directory
files = dir('SUBJ_LEMON_02');
files=files(5:end);
addpath(files(1).folder)

split_percentage=0.2;
%create struct
ttest_info=struct('pvals',0,'significant',0,'ptrain',0,'corrtrain',0);

%create covariate labels
for i=0:length(node_labels)-1
        for j=1:length(labels_hrf)
             COVNAMES(floor(i*9+j))=strcat(node_labels(i+1),'_',labels_hrf(j));
        end
end

%% save univariate selection across 1000 iterations

for ss=1:1000
    
    ss
    
    %set random seed
    rng(ss,'Threefry');
   
    [rsn_hrfy,rsn_hrfo,Y_LR_y,Y_LR_o,S_final,partition,isSignificant,adjusted_pvals]=...
        univariate_selection(files,node_labels,labels_hrf,ss,split_percentage,node_bounds);
    
    ttest_info(ss).pvals=adjusted_pvals;
    ttest_info(ss).significant=isSignificant;
    
    %% compute and save data correlation with output (for data representativeness)

    %define X and Y
    yo=reshape(rsn_hrfy,666,size(Y_LR_y,2));
    old=reshape(rsn_hrfo,666,size(Y_LR_o,2));
    X = [yo,old]';     %columns of covariates
    Y = [ones(size(Y_LR_y,2),1); zeros(size(Y_LR_o,2),1)]; %columns of outputs
    Y_LR=[Y_LR_y,Y_LR_o]';
    
    %stratified Training-Test split
    Train_idx = training(partition);
    Test_idx = test(partition);
    
    % split the input
    XTrain = X(Train_idx,selec);
    XTest = X(Test_idx,selec);
    
  
    % split the output
    YTrain = Y(Train_idx,:);
    YTest = Y(Test_idx,:);
    Y_LRTrain=Y_LR(Train_idx,:);
    Y_LRTest=Y_LR(Test_idx,:);
    
    %correlation with output
    [R_train,pval_train]=corr(X(Train_idx,:),Y_LRTrain);
    
    ttest_info(ss).ptrain=pval_train;
    ttest_info(ss).corrtrain=R_train;

  


end


