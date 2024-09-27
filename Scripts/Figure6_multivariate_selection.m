clear all
clc
close all

load node_lab_conc.mat
load group_scatter.mat

%load correlation with all dataset
load CORRg.mat

labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];
labels_hrf_plot=[{'A peak'},{'A trough'},{'t peak'},{'t trough'},{'Rise slope'},{'Fall slope'},...
    {'AUC'},{'FWHM'},{'Peak to trough'}];

base_path=pwd;
files=dir(fullfile(base_path,'FINAL_RESULT_folder')); %data obtainedby running main_multicollinearity_multivariate_selection.m
addpath(genpath(files(1).folder))
files=files(4:end);

%create data mat
data=[];
for ii=1:length(files)
    load(files(ii).name);
    data=[data,FINAL_RESULT];
    clear FINAL_RESULT
end

%create covariate labels
for ii=0:length(node_labels)-1
        for j=1:length(labels_hrf)
             CovNames(floor(ii*9+j))=strcat(node_labels(ii+1),'_',labels_hrf(j));
        end
end

%initialize empty mat variables
selec_sum=zeros(length(data),length(CovNames));
selec_beta=zeros(length(data),length(CovNames));
selec_zscore=zeros(length(data),length(CovNames));

%% across 1000 iterations

for ii=1:length(data)

    complexity(ii)=length(data(ii).Covariates);
    selec_name=data(ii).Coeffs.Properties.RowNames(2:end);

    if isempty(selec_name)
        selec=zeros(length(CovNames));
        continue;
    end

    selec=double(matches(CovNames,selec_name));
    selec_sum(ii,:)=selec;
    intercept(ii)=data(ii).Coeffs.Estimate(1);
    intercept_zscore(ii)=data(ii).Coeffs.pValue(1);
    selec_beta(ii,find(selec))=data(ii).Coeffs.Estimate(2:end);
    selec_zscore(ii,find(selec))=data(ii).Coeffs.pValue(2:end);
    
    %compute distance between correlation vectors
    dist_corr(ii)=pdist([R,data(ii).corrtrain]','euclidean'); 
    map(ii,complexity(ii)-1)=dist_corr(ii);
end

%select only good models
AUC_iter=[data.AUC];
AUC_iter_beta=repmat(AUC_iter',1,length(CovNames));
complexity_beta=repmat(complexity',1,length(CovNames));
intercept_clean=intercept.*(intercept_zscore<0.05 & intercept_zscore>0).*(AUC_iter>0.5).*(complexity>=4 & complexity<=6);
selec_beta_clean=selec_beta.*(selec_zscore<0.05 & selec_zscore>0).*(AUC_iter_beta>0.5).*(complexity_beta>=4 & complexity_beta<=6);

%selection frequency
SUM_cov_post=(sum((selec_beta_clean~=0))/1000)*100;

for ll=1:length(labels_hrf)

    HRF_perc_selection2_post(:,ll)=SUM_cov_post(ll:9:end);

end

%reorder subcorticals and exclude cerebellum
HRF_perc_selection2_post_LR=HRF_perc_selection2_post([1:66,68:73],:);
HRF_perc_selection2_post_LR_ord=HRF_perc_selection2_post_LR;
HRF_perc_selection2_post_LR_ord([63:end],:)=HRF_perc_selection2_post_LR([63 63+5 64 64+5 65 65+5 66 66+5 67 67+5],:);

node_labels_S=load('node_labels_Schaefer.mat').node_labels;
node_lab_LR=node_labels_S([1:66,68:73]);
node_lab_LR_ord=node_lab_LR;
node_lab_LR_ord([63:end])=node_lab_LR([63 63+5 64 64+5 65 65+5 66 66+5 67 67+5]);


%% frequency plot

%RSN color coding
group_col = [0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
    0.466666666666667   0.674509803921569   0.188235294117647;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
    0.929411764705882   0.694117647058824   0.125490196078431;0.850980392156863   0.325490196078431   0.098039215686275;...
         0    1.0000    0.5000
    1.0000         0    0.5000
    0.7500    0.7500         0
    1.0000         0         0
         0         0    0.5000];


group_nets=[1:7,1:7,8:12];

marker_list={'o','+','x','s','^','diamond','v','pentagram','>'};
color_list=gray(9);
bounds_x=[0 2 5 9 16 19 23 32 35 38 43 47 49 55 62 64 66 68 70]+0.5;
bounds_y=[-0.5 -0.5 110 110];

figure
for ll=1:length(labels_hrf)
    stem(HRF_perc_selection2_post_LR_ord(:,ll),'Marker',marker_list(ll),'MarkerEdgeColor','k','MarkerFaceColor',color_list(ll,:),...
        'LineStyle',':','Color','k','MarkerSize',7)
    hold on
end
legend(labels_hrf_plot,'Interpreter','none')

hold on
for nn=1:length(bounds_x)-1
    patch([bounds_x(nn) bounds_x(nn+1) bounds_x(nn+1) bounds_x(nn)],bounds_y,group_col(group_nets(nn),:),...
        'FaceAlpha',0.2)
    hold on
end

nn=length(bounds_x);
patch([bounds_x(nn) length(node_lab_LR_ord)+0.5 length(node_lab_LR_ord)+0.5 bounds_x(nn)],bounds_y,group_col(group_nets(nn),:),'FaceAlpha',0.1)
    
xticks(1:72)
set(gca,'XTickLabel',node_lab_LR_ord,'XTickLabelRotation',90)
set(gca,'TickLength',[0 0])
ylim([-0.5 105])
xlim([0.5 72.5])
ylabel('Selection frequency (%)')

%% frequency descending plot

SUM_cov_univar_sel=load('perc_selection1.mat').perc_selection1;
[feat_ordered,idx_ordered]=sort(SUM_cov_post,'descend'); 
selec_final_feat=sort(idx_top(1:5));

for ii=0:length(node_labels_S)-1
        for j=1:length(labels_hrf)
             CovNames_new(floor(ii*9+j))=strcat(node_labels_S(ii+1),'_',labels_hrf(j));
             CovNames_struct(floor(ii*9+j)).Cov=CovNames_new(floor(ii*9+j));
             CovNames_struct(floor(ii*9+j)).ROI=node_labels_S(ii+1);
             CovNames_struct(floor(ii*9+j)).HRF_param=labels_hrf(j);
        end
end

figure
%select top 10 features
for jj=1:10
    
    Covariate_name=CovNames_struct(idx_ordered(jj)).Cov;
    Covariate_ROI=CovNames_struct(idx_ordered(jj)).ROI;
    Covariate_HRF_param=CovNames_struct(idx_ordered(jj)).HRF_param;
    idx_HRF_col=find(matches(labels_hrf,Covariate_HRF_param));
    idx_ROI_col=find(matches(node_lab_LR,Covariate_ROI));

    
    stem(jj-0.1,log(SUM_cov_post(idx_ordered(jj))),'Marker',marker_list(idx_HRF_col),'MarkerFaceColor',color_list(idx_HRF_col,:),...
        'MarkerSize',9,'MarkerEdgeColor','k','LineStyle',':','Color','k')
    x_plot1(jj)=jj-0.1;
    y_plot1(jj)=log(SUM_cov_post(idx_ordered(jj)));
    hold on

    stem(jj+0.1,log(SUM_cov_univar_sel(idx_ordered(jj))),'Marker',marker_list(idx_HRF_col),'MarkerFaceColor',color_list(idx_HRF_col,:),...
        'MarkerSize',9,'MarkerEdgeColor','k','LineStyle',':','Color','k')
    xlim([0 11.5])
    tick_col=[0 0 0];
    ticklabel_col{jj}=['\color[rgb]{' num2str(tick_col(1)) ',' num2str(tick_col(2)) ',' num2str(tick_col(3)) '}' Covariate_ROI{1}];
  
    hold on


end

set(gca,'XTick',1:length(idx_ordered), 'XTickLabel', ticklabel_col,'XTickLabelRotation',90);
set(gca,'TickLength',[0 0])
ylabel('Selection frequency % (log scale)')
hold on
plot(x_plot1,y_plot1,'b:','LineWidth',2)
hold on
ylimits=get(gca,'YLim');
patch([0 5.5 5.5 0],[ylimits(1) ylimits(1) ylimits(2) ylimits(2)],'b','FaceAlpha',0.1)
hold off

