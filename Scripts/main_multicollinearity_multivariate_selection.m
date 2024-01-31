clear all
clc
close all

%% load data and initialize variables
load node_lab_conc.mat
labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];
node_bounds=1:74;

files = dir('...');
files=files(3:end);
addpath(files(1).folder)
split_percentage=0.2;
FINAL_RESULT=struct();

%covariate names
for i=0:length(node_labels)-1
        for j=1:length(labels_hrf)
             COVNAMES(floor(i*9+j))=strcat(node_labels(i+1),'_',labels_hrf(j));
        end
end

%struct to save results from multicollinearity removal
corr_step=struct();

%% fit 1000 model realizations

parfor ss=1:1000 

    %set seed for reproducibility
    seed = ss
    rng(seed,'Threefry');
    perm=0;
    
    %% univariate selection
    [rsn_hrfy,rsn_hrfo,Y_LR_y,Y_LR_o,S_final,partition,isSignificant,adjusted_pvals]=univariate_selection(files,node_labels,labels_hrf,seed,split_percentage,node_bounds);
    CovNames=COVNAMES;
    selec_name=S_final;
    selec=find(matches(CovNames,selec_name));
    CovNames=CovNames(:,selec);

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
    
    % compute correlation with output (for data representativeness)
    [R_train,pval_train]=corr(X(Train_idx,:),Y_LRTrain);


    %% multicollinearity removal step
    
    S=1;

    %set condition number
    [sValue,condIdx,VarDecomp] = collintest(XTrain);
    CN=condIdx(end);

    %remove redundant features
    while S~=0
        thr_corr=0.85;
        CorrMat = corr(XTrain,'rows','pairwise');
        mat_thr=CorrMat.*(abs(CorrMat)>thr_corr)-diag(diag(eye(length(CovNames))));
        id_corr=sum(mat_thr);
        S=sum(sum(id_corr));
        if S==0
            break;
        end
   
        [removed_i,removed_j] =(find(mat_thr~=0,1,'first'));

        %retain the most correlated feature
        if abs(corr(XTrain(:,removed_i),Y_LRTrain))>abs(corr(XTrain(:,removed_j),Y_LRTrain))
            cov_index = not(matches(CovNames,CovNames(removed_j)));
        else
            cov_index = not(matches(CovNames,CovNames(removed_i)));
        end
    
        XTrain = XTrain(:,cov_index);
        %condition number
        [sValue,condIdx,VarDecomp] = collintest(XTrain);
        CN_iter=condIdx(end);
        CN=[CN,CN_iter];

        XTest = XTest(:,cov_index);
        CovNames = CovNames(:,cov_index);
        
    end
    
    corr_step(ss).Cov=CovNames;
    corr_step(ss).CN=CN;

    %% Bootstrap sampling 
    
    
    k=100; % number of iterations

    %save original sets
    XTrain_orig=XTrain;
    YTrain_orig=YTrain;
    YTest_orig=YTest;
    XTest_orig=XTest;
    
    weights_train_1=sum(YTrain_orig);
    weights_train_0=length(YTrain_orig)-weights_train_1;
    prob_0=weights_train_0/length(YTrain_orig);
    prob_1=1-prob_0;
   
    %Initialize empty struct to be used through CV
    CV_DATA=struct();
    
    
    for f=1:k 
            
  
        check_Train=0;
        check_Test=0;

        %check balance of classes in the train and test subsets (same
        %percentage of old and young subjects of the original dataset)
        while check_Train>prob_1+0.02 || check_Train<prob_1-0.02 || check_Test>prob_1+0.02 || check_Test<prob_1-0.02
        
            [w,Train_idx_cv]=datasample(YTrain_orig,length(YTrain_orig));
            Test_idx_cv=setdiff([1:length(YTrain_orig)],Train_idx_cv);
            XTrain=XTrain_orig(Train_idx_cv,:);
            YTrain=YTrain_orig(Train_idx_cv);
            XTest=XTrain_orig(Test_idx_cv,:);
            YTest=YTrain_orig(Test_idx_cv);
            check_Train=sum(YTrain)/length(Train_idx_cv)
            check_Test=(sum(YTest)/length(Test_idx_cv))
            
        end  
           
           
        %variable normalization
        XTrain = (XTrain-mean(XTrain))./std(XTrain);
        XTest  = (XTest-mean(XTrain))./std(XTrain);
        
            
           
        %% Stepwise backward selection
    
    
        varNames = matlab.lang.makeValidName([CovNames,{'Output'}]);
        warning off
        step_model = stepwiseglm(XTrain, YTrain, 'linear', 'Distribution', 'binomial',...
            'link', 'logit', 'upper', 'linear', 'criterion', 'bic', ...
            'DispersionFlag',true,'VarNames',varNames)%,'Lower','linear');
        warning on
    
        %save selected features in CV_DATA
        CV_DATA(f).features_S=step_model.PredictorNames;
        
     
    end


    %% Check for features selected across bootstrap samples
     
    %check how many times features have been selected 
    Cov_percent_S=zeros(size(CovNames));
    for f=1:k
        
        for ll=1:length(CovNames)
            if sum(matches(CV_DATA(f).features_S,CovNames(ll)))==1
               Cov_percent_S(ll)=Cov_percent_S(ll)+1;
            end
        end
    end
    
    fprintf('\n\tVARIABLE\tNUMBER OF FOLDS(STEPWISE)\n')
    for i=1:length(CovNames)
       fprintf('%12s\t%+1f\n',string(CovNames(i)),Cov_percent_S(i)) 
    end    
    
    
    %% TRAINING of final model (for iteration ss)
    
    XTrain=XTrain_orig;
    YTrain=YTrain_orig;
    XTest=XTest_orig;
    YTest=YTest_orig;
    
    
    %variable normalization
    XTrain = (XTrain-mean(XTrain))./std(XTrain);
    XTest  = (XTest-mean(XTrain))./std(XTrain);
    
    
    %% STEPWISE model with features selected >=50%

    thr_sel=floor((k/100)*50);
    varNames_S = matlab.lang.makeValidName([CovNames(Cov_percent_S>=thr_sel),{'Output'}]);
    step_model=fitglm(XTrain(:,Cov_percent_S>=thr_sel), YTrain, 'linear', 'Distribution', 'binomial',...
            'link', 'logit','DispersionFlag',true,'VarNames',varNames_S);
   
    disp(step_model)
    
    
    %% TEST of final model
    
    %indices for feature selection on Test set
    valid_idx_S=matches(CovNames,step_model.PredictorNames);
    XTest_S=XTest(:,valid_idx_S);
    b_S=step_model.Coefficients.Estimate; 
    p_S=glmval(b_S,XTest_S,'logit');

    %Likelihood ratio-R^2 statistics
    R2_lr_all=step_model.Rsquared.Deviance;
    %Cox&Snell-R^2 statistics
    R2_cs_all=step_model.Rsquared.AdjGeneralized;
    %McFadden's-R^2 statistics
    R2_mcf_all=step_model.Rsquared.LLR;
    %C-index
    C_I_all=c_index(p_S,YTest);
    
    
   
    fprintf('SUMMARY STATISTICS FOR STEPWISE MODEL \n')
    disp(['Likelihood ratio-R2 statistic: ',num2str(R2_lr_all)])
    disp(['Cox&Snell-R2 statistic: ',num2str(R2_cs_all)])
    disp(['McFadden-R2 statistic: ',num2str(R2_mcf_all)])
    disp(['Concordance Index: ',num2str(C_I_all)])
    %Deviance
    disp(['Deviance: ',num2str(step_model.Deviance)])
    %AIC and BIC
    disp(['AIC: ',num2str(step_model.ModelCriterion.AIC)])
    disp(['BIC: ',num2str(step_model.ModelCriterion.BIC)])
    
    

    FINAL_RESULT(ss).Covariates=step_model.CoefficientNames;
    FINAL_RESULT(ss).Coeffs=step_model.Coefficients;
    FINAL_RESULT(ss).AUC=C_I_all;
    FINAL_RESULT(ss).AIC=step_model.ModelCriterion.AIC;
    FINAL_RESULT(ss).R2_CS=R2_cs_all;
    FINAL_RESULT(ss).corrtrain=R_train;
    FINAL_RESULT(ss).ptrain=pval_train;


end

