clear all
clc
close all

addpath(genpath('ENIGMA')) %load ENIGMA toolbox
addpth('Functions')
load node_lab_conc.mat
load group_scatter.mat

labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},...
    {'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];

node_bounds=1:74; %cortical+subcortical ROIs
TR=1.4; %s

%load data directory
files = dir('SUBJ_LEMON_02');
files=files(5:end);
addpath(files(1).folder)

%create covariate labels
for i=0:length(node_labels)-1
        for jj=1:length(labels_hrf)
             COVNAMES(floor(i*9+jj))=strcat(node_labels(i+1),'_',labels_hrf(jj));
        end
end

%% young vs old HRF partition

age_young=30;
age_old=60;

%peforms HRF partition into young vs old groups
[rsn_hrfy,rsn_hrfo,HRFmatrix_y,HRFmatrix_o,Y_LR_y,Y_LR_o] = HRF_partition(files,age_young,age_old);
    

%% only cortical ROIs (first 62 ROIs)

HRFmatrix_o=HRFmatrix_o(:,1:62,:);
HRFmatrix_y=HRFmatrix_y(:,1:62,:);
group_num=group_num(1:62);
node_labels=node_labels(1:62);

hrf_young=reshape(HRFmatrix_y,size(HRFmatrix_y,1),length(node_labels)*size(HRFmatrix_y,3));
hrf_old=reshape(HRFmatrix_o,size(HRFmatrix_o,1),length(node_labels)*size(HRFmatrix_o,3));

%PCA dimensionality reduction
DATA=([hrf_young,hrf_old]'); 
HRF_length=size(HRFmatrix_y,1);

[tot_perc,pc1_young,pc2_young,pc3_young,pc1_old,pc2_old,pc3_old] = PCA_analysis(DATA,TR,HRF_length,hrf_young,hrf_old);


%RSN color coding
group_col = [0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
    0.466666666666667   0.674509803921569   0.188235294117647;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
    0.929411764705882   0.694117647058824   0.125490196078431;0.850980392156863   0.325490196078431   0.098039215686275];

group_rep_young=repmat(group_num,1,size(HRFmatrix_y,3));
group_rep_old=repmat(group_num,1,size(HRFmatrix_o,3));
group_rep_nodes_young=repmat(1:62,1,size(HRFmatrix_y,3));
group_rep_nodes_old=repmat(1:62,1,size(HRFmatrix_o,3));

%select example ROI
node1=53;
find_node_young1=group_rep_nodes_young==node1;
find_node_old1=group_rep_nodes_old==node1;


figure
plot(cumsum(explained),'.-','MarkerSize',12,'MarkerFaceColor','b','Color','k')
xlabel('Number of PCs')
ylabel('Percentage variance explained')
title('PCA cumulative scree plot')

figure
gscatter3(pc1_young,pc2_young,pc3_young,group_rep_young,group_col,{'s','s','s','s','s','s','s'},5)
hold on
scatter3(pc1_young(find_node_young1),pc2_young(find_node_young1),pc3_young(find_node_young1),'s','MarkerEdgeColor','k','MarkerFaceColor',group_col(group_num(node1),:),'LineWidth',0.7)
hold on
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
xlimits=get(gca,'XLim');
ylimits=get(gca,'YLim');
zlimits=get(gca,'ZLim');

figure
gscatter3(pc1_old,pc2_old,pc3_old,group_rep_old,group_col,{'s','s','s','s','s','s','s'},5)
hold on
scatter3(pc1_old(find_node_old1),pc2_old(find_node_old1),pc3_old(find_node_old1),'s','MarkerEdgeColor','k','MarkerFaceColor',group_col(group_num(node1),:),'LineWidth',0.7)
hold on
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
xlim(xlimits)
ylim(ylimits)
zlim(zlimits)


% display HRFs
F=0:TR:size(HRFmatrix_y,1)*TR-TR;
figure()
    
find_nodes_young=group_rep_nodes_young==node1;
stdshade(hrf_young(:,find_nodes_young)',0.2,[0.074509803921569   0.623529411764706   1.000000000000000],F)
xlabel('Time (s)')
ylimits=get(gca,'YLim');
hold on
find_nodes_old=group_rep_nodes_old==node1;
stdshade(hrf_old(:,find_nodes_old)',0.2,[0.850980392156863   0.325490196078431   0.098039215686275],F)


% Map parcellated data to the surface
corr_map_conte69 = parcel_to_surface(group_num,'schaefer_62_conte69');
% Project the results on the surface brain
g= figure,plot_cortical_Schaefer(corr_map_conte69, 'surface_name', 'conte69')



%% only subcortical ROIs

%exclude cerebellum from analysis
HRFmatrix_o=HRFmatrix_o(:,[63:66,68:73],:);
HRFmatrix_y=HRFmatrix_y(:,[63:66,68:73],:);
group=group([63:66,68:73]);

group_num=[1 2 3 4 5 1 2 3 4 5];
node_labels=node_labels([63:66,68:73]);

hrf_young=reshape(HRFmatrix_y,size(HRFmatrix_y,1),length(node_labels)*size(HRFmatrix_y,3));
hrf_old=reshape(HRFmatrix_o,size(HRFmatrix_o,1),length(node_labels)*size(HRFmatrix_o,3));

%PCA dimensionality reduction
DATA=([hrf_young,hrf_old]'); 

[tot_perc,pc1_young,pc2_young,pc3_young,pc1_old,pc2_old,pc3_old] = PCA_analysis(DATA,TR,HRF_length,hrf_young,hrf_old);


%RSNs color coding
reorder=[2 5 4 3 1];
group_col_sub=[1 0 0.5; 0 0 0.5;...
    1 0 0;0.75 0.75 0;0 1 0.5];

group_col_sub(reorder,:)=group_col_sub;

%select example ROI
node1=1;
group_rep_young=repmat(group_num,1,size(HRFmatrix_y,3));
group_rep_old=repmat(group_num,1,size(HRFmatrix_o,3));
group_rep_nodes_young=repmat(1:10,1,size(HRFmatrix_y,3));
group_rep_nodes_old=repmat(1:10,1,size(HRFmatrix_o,3));

find_node_young1=group_rep_nodes_young==node1;
find_node_old1=group_rep_nodes_old==node1;


figure
plot(cumsum(explained),'.-','MarkerSize',12,'MarkerFaceColor','b','Color','k')
xlabel('Number of PCs')
ylabel('Percentage variance explained')
title('PCA cumulative scree plot')

figure
gscatter3(pc1_young,pc2_young,pc3_young,group_rep_young,group_col_sub,{'s','s','s','s','s'},5)
hold on
scatter3(pc1_young(find_node_young1),pc2_young(find_node_young1),pc3_young(find_node_young1),'s','MarkerEdgeColor','k','MarkerFaceColor',group_col_sub(group_num(node1),:),'LineWidth',0.7)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
xlimits=get(gca,'XLim');
ylimits=get(gca,'YLim');
zlimits=get(gca,'ZLim');

figure
gscatter3(pc1_old,pc2_old,pc3_old,group_rep_old,group_col_sub,{'s','s','s','s','s'},5)
hold on
scatter3(pc1_old(find_node_old1),pc2_old(find_node_old1),pc3_old(find_node_old1),'s','MarkerEdgeColor','k','MarkerFaceColor',group_col_sub(group_num(node1),:),'LineWidth',0.7)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
xlim(xlimits)
ylim(ylimits)
zlim(zlimits)


%display HRFs
F=0:TR:size(HRFmatrix_y,1)*TR-TR;

figure
find_nodes_young=group_rep_nodes_young==node1;
stdshade(hrf_young(:,find_nodes_young)',0.2,[0.074509803921569   0.623529411764706   1.000000000000000],F)
xlabel('Time (s)')
ylimits=get(gca,'YLim');
hold on
find_nodes_old=group_rep_nodes_old==node1;
stdshade(hrf_old(:,find_nodes_old)',0.2,[0.850980392156863   0.325490196078431   0.098039215686275],F)

%surface plot
subcortical_idx=group_num;
figure
g_sub=plot_subcortical_AAL3(subcortical_idx,'ventricles','False')


