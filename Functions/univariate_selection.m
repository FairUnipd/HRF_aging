function [rsn_hrfy,rsn_hrfo,Y_LR_y,Y_LR_o,S_final,partition,isSignificant,adjusted_pvals]=univariate_selection(files,node_labels,labels_hrf,seed,split_percentage,node_bounds)
 
    %define counters
    county=1;
    counto=1;
    

    age_young=30; %young cut-off
    age_old=60;   %old cut-off
    
    num_HRF_feat=9;

    %define covariate names
    for ii=0:length(node_labels)-1
        for jj=1:length(labels_hrf)
             CovNames(floor(ii*num_HRF_feat+jj))=strcat(node_labels(ii+1),'_',labels_hrf(jj));
        end
    end
    
    %peforms HRF partition into young vs old groups
    [rsn_hrfy,rsn_hrfo,~,~,Y_LR_y,Y_LR_o] = HRF_partition(files,age_young,age_old);
    
    
    % apply hold-out dataset partition
    rng(seed)
    Y = [ones(size(Y_LR_y,2),1); zeros(size(Y_LR_o,2),1)]; %columns of outputs
    group = Y;
    partition = cvpartition(group,'Holdout',split_percentage,'Stratify',true);
    

    %statistical one-sample t-test

    for ii=1:size(rsn_hrfo,1) %for each parameter

        for jj=1:length(node_bounds) %for each ROI
            rsn_hrfy_ij=rsn_hrfy(ii,jj,:);
            rsn_hrfo_ij=rsn_hrfo(ii,jj,:);
            

            [hy,pval_rsn(ii,jj)] = ttest2(rsn_hrfy_ij,rsn_hrfo_ij,'Vartype','unequal');
            if isnan(pval_rsn(ii,jj))
                pval_rsn(ii,jj)=1;
            end
        end
    end

    % Bonferroni-Holm correction
    pvals=pval_rsn(:);
    disp((sum(isnan(pvals))))
    
    [isSignificant,adjusted_pvals,alpha]= bonferroni_holm(pvals,0.05);
    
    
    isSignificant_rsn=reshape(isSignificant(1:end),num_HRF_feat,length(node_bounds));
    cs=1;
    S=struct([]); %save indices of significant features
    for jj=1:length(node_bounds)
        for ii=1:size(rsn_hrfo,1)
            if isSignificant_rsn(ii,jj)==1
                
                S(cs).young=squeeze(rsn_hrfy(ii,jj,:));
                S(cs).old=squeeze(rsn_hrfo(ii,jj,:));
                S(cs).label=ROI_labels(jj);
                S(cs).param=labels_hrf(ii);
                cs=cs+1;
               
            end
        
        end
    end
   

    for ii=1:length(S)

        if isempty(S(1)) %if no covariate is selected

            selec_name={};
            break;
        end
        selec_name(ii)=strcat(S(ii).label,'_',S(ii).param);
    end
    
    %discard useless covariates
    selec=find(matches(CovNames,selec_name));
    idx_k=zeros(size(CovNames));
    idx_k(selec)=1;
   
    %select significant covariates
    S_final=CovNames(logical(idx_k))';

end
