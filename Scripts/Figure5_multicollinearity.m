clear all
clc
close all

addpath('Data')
addpath('Functions')

%load data
load('AFTER_CORR_STEP_DATA.mat')
load node_lab_conc.mat
load group_scatter.mat %RSN labels

%HRF features labels
labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];
labels_hrf_plot=[{'A peak'},{'A trough'},{'t peak'},{'t trough'},{'Rise slope'},{'Fall slope'},{'AUC'},{'FWHM'},{'Peak to trough'}];

%% identify selected features

%initial covariates
for ii=0:length(node_labels)-1
        for j=1:length(labels_hrf)
             CovNames(floor(ii*9+j))=strcat(node_labels(ii+1),'_',labels_hrf(j));
        end
end

%analyse iterative results
for ii=1:length(corr_step)

    selec_name=corr_step(ii).Cov;

    if isempty(selec_name)
        selec=zeros(length(CovNames));
        continue;
    end

    selec=double(matches(CovNames,selec_name));
    feat_selection2(ii,:)=selec;
    
end

perc_selection2=(sum(feat_selection2)/1000)*100;

for ll=1:length(labels_hrf)

    HRF_perc_selection2(:,ll)=perc_selection2(ll:9:end);

end

%reorder subcorticals
node_labels=load('node_labels_Schaefer.mat').node_labels;
HRF_perc_selection2_LR=HRF_perc_selection2([1:66,68:73],:);
HRF_perc_selection2_LR_ord=HRF_perc_selection2_LR;
HRF_perc_selection2_LR_ord([63:end],:)=HRF_perc_selection2_LR([63 63+5 64 64+5 65 65+5 66 66+5 67 67+5],:);

node_lab_LR=node_labels([1:66,68:73]);
node_lab_LR_ord=node_lab_LR;
node_lab_LR_ord([63:end])=node_lab_LR([63 63+5 64 64+5 65 65+5 66 66+5 67 67+5]);


%% Panel A: feature selection figure

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
    stem(HRF_perc_selection2_LR_ord(:,ll),'Marker',marker_list(ll),'MarkerEdgeColor','k','MarkerFaceColor',color_list(ll,:),...
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
ylim([-1 105])
xlim([0.5 72.5])
ylabel('Selection frequency (%)')
hold off

%% Panel B: Condition number plot
count=0;

for ii=1:length(corr_step)

    num_iter=length(corr_step(ii).CN);
    CN_iter=corr_step(ii).CN;
    
    %interpolation for better visualization
    x_interp=[1:1:10];
    if num_iter<=1
        continue
    end
    count=count+1;
    tmp=griddedInterpolant(1:num_iter,CN_iter,'linear','linear');
    CN_interp(:,count)=tmp(x_interp);
    
end

%box plot graph
figure
hold on
colors=winter(max(x_interp));
N=size(CN_interp,2);
group_box=[ones(1,N)';2*ones(1,N)';3*ones(1,N)';4*ones(1,N)';5*ones(1,N)';6*ones(1,N)';7*ones(1,N)';8*ones(1,N)';9*ones(1,N)';10*ones(1,N)'];
DATA_box=[];
for rr=1:max(x_interp)
    DATA_box=[DATA_box;CN_interp(rr,:)'];
end
pos=[1:max(x_interp)];
symbol={'o','o','o','o','o','o','o','o','o','o'};
plot_box_scatter(DATA_box,group_box,pos,colors,symbol)
hold on
ylim([-700 800])
xlim([0.5 10.5])
xlabel('iterations')
ylabel('CN')