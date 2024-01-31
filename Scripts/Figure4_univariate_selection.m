clear all
clc
close all

load('info_ttest.mat')
load node_labels_Schaefer.mat
load group_scatter.mat
labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];

labels_hrf_4plot=[{'A peak'},{'A trough'},{'t peak'},{'t trough'},{'Rise slope'},{'Fall slope'},{'AUC'},{'FWHM'},{'Peak to trough'}];

%% compute percentage of selection
feat_selection1=[ttest_info.significant];
perc_selection1=(sum(feat_selection1')/1000)*100;

for ll=1:length(labels_hrf)

    HRF_perc_selection1(:,ll)=perc_selection1(ll:9:end);

end

%% reorder subcortical and exclude cerebellum (for visualization)
HRF_perc_selection1_LR=HRF_perc_selection1([1:66,68:73],:);
HRF_perc_selection1_LR_ord=HRF_perc_selection1_LR;
HRF_perc_selection1_LR_ord([63:end],:)=HRF_perc_selection1_LR([63 63+5 64 64+5 65 65+5 66 66+5 67 67+5],:);

node_lab_LR=node_labels([1:66,68:73]);
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
    stem(HRF_perc_selection1_LR_ord(:,ll),'Marker',marker_list(ll),'MarkerEdgeColor','k','MarkerFaceColor',color_list(ll,:),...
        'LineStyle',':','Color','k','MarkerSize',7)
    hold on
end
legend(labels_hrf_4plot,'Interpreter','none')

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


%% stem plots
load HRF_median_maps_no_cer_not_reordered.mat

colors_maps=color_list;

figure
for ll=1:length(labels_hrf)

    subplot(9,1,ll)
    map_young_cort=HRF_maps(ll).cortical_map_young;
    map_old_cort=HRF_maps(ll).cortical_map_old;
    map_young_subcort=HRF_maps(ll).subcortical_map_young;
    map_old_subcort=HRF_maps(ll).subcortical_map_old;
    map_young=[map_young_cort,map_young_subcort];
    map_old=[map_old_cort,map_old_subcort];
    map_diff=map_young-map_old;
    map_diff=map_diff.*(HRF_perc_selection1_LR(:,ll)>0)';
    vec=1:72;
    vec=vec.*(HRF_perc_selection1_LR(:,ll)>0)';
    label=HRF_maps(ll).label;
    stem(vec,map_diff,'Marker',marker_list(ll),'MarkerEdgeColor','k','MarkerFaceColor',color_list(ll,:),...
        'LineStyle',':','Color','k','MarkerSize',7)
    hold on
    
    if ll==length(labels_hrf)
        xticks(1:72)
        set(gca,'XTickLabel',node_lab_LR,'XTickLabelRotation',90)
        set(gca,'TickLength',[0 0])
    end
    ylim([-max(abs(map_diff))-0.1*max(abs(map_diff)) max(abs(map_diff))+0.1*max(abs(map_diff))])
    xlim([1 72])

    patch([1 72 72 1],[0 0 max(abs(map_diff))+0.1*max(abs(map_diff)) max(abs(map_diff))+0.1*max(abs(map_diff))],...
        [0.074509803921569   0.623529411764706   1.000000000000000],'FaceAlpha',0.2)
    hold on
    patch([1 72 72 1],[0 0 -max(abs(map_diff))-0.1*max(abs(map_diff)) -max(abs(map_diff))-0.1*max(abs(map_diff))],...
        [0.850980392156863   0.325490196078431   0.098039215686275],'FaceAlpha',0.2)
 
end