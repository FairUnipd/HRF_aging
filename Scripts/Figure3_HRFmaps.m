clear all
clc
close all

%load ENIGMA toolbox
addpath(genpath('ENIGMA'))
load node_lab_conc.mat
load group_scatter.mat
labels_hrf=[{'A_peak'},{'A_trough'},{'t_peak'},{'t_trough'},{'Rise_slope'},{'Fall_slope'},{'AUC'},{'FWHM'},{'Peak_to_trough'}];
node_bounds=1:74;

%load data directory
files = dir('SUBJ_LEMON_02');
files=files(5:end);
addpath(files(1).folder)

%covariate labels
for i=0:length(node_labels)-1
        for j=1:length(labels_hrf)
             COVNAMES(floor(i*9+j))=strcat(node_labels(i+1),'_',labels_hrf(j));
        end
end

%% compute HRFs curves

age_young=30;
age_old=60;

%peforms HRF partition into young vs old groups
[params_y,params_o,HRFmatrix_y,HRFmatrix_o,Y_LR_y,Y_LR_o] = HRF_partition(files,age_young,age_old);
  

%% median maps of HRF parameters

%create HRF parametric median maps
HRF_maps=struct();
for pp=1:length(labels_hrf)
    HRF_maps(pp).label=labels_hrf(pp);
    HRF_maps(pp).cortical_map_young=median(squeeze(params_y(pp,1:62,:))');
    HRF_maps(pp).subcortical_map_young=median(squeeze(params_y(pp,[63:66,68:73],:))');
    HRF_maps(pp).cortical_map_old=median(squeeze(params_o(pp,1:62,:))');
    HRF_maps(pp).subcortical_map_old=median(squeeze(params_o(pp,[63:66,68:73],:))');
end

%% young subjects

for idx_feat=1:length(labels_hrf)

    % Map parcellated data to the surface
    HRF_map_cortical = parcel_to_surface(HRF_maps(idx_feat).cortical_map_young,'schaefer_62_conte69');

    % Project the results on the surface brain
    if idx_feat==3
        range_cort=[1 7];
        range_subcort=[1 7];
    else
        range_cort=[min([HRF_maps(idx_feat).cortical_map_young HRF_maps(idx_feat).cortical_map_old]) max([HRF_maps(idx_feat).cortical_map_young HRF_maps(idx_feat).cortical_map_old])];
        range_subcort=[min([HRF_maps(idx_feat).subcortical_map_young HRF_maps(idx_feat).subcortical_map_old]) max([HRF_maps(idx_feat).subcortical_map_young HRF_maps(idx_feat).subcortical_map_old])+10e-10];
    end
    
    %corticals 
    g= figure,plot_cortical(HRF_map_cortical, 'surface_name', 'conte69','cmap','RdBu_r','color_range',range_cort)
    sgtitle([labels_hrf(idx_feat) ' young'],'Interpreter','none')

    %subcorticals
    g_sub=figure, plot_subcortical_AAL3(HRF_maps(idx_feat).subcortical_map_young,'ventricles','False', 'cmap', 'RdBu_r','color_range',range_subcort)
    sgtitle([labels_hrf(idx_feat) ' young'],'Interpreter','none')

end

%% old subjects

for idx_feat=1:length(labels_hrf)

    % Map parcellated data to the surface
    HRF_map_cortical = parcel_to_surface(HRF_maps(idx_feat).cortical_map_old,'schaefer_62_conte69');

    % Project the results on the surface brain
    if idx_feat==3
        range_cort=[1 7];
        range_subcort=[1 7];
    else
        range_cort=[min([HRF_maps(idx_feat).cortical_map_young HRF_maps(idx_feat).cortical_map_old]) max([HRF_maps(idx_feat).cortical_map_young HRF_maps(idx_feat).cortical_map_old])];
        range_subcort=[min([HRF_maps(idx_feat).subcortical_map_young HRF_maps(idx_feat).subcortical_map_old]) max([HRF_maps(idx_feat).subcortical_map_young HRF_maps(idx_feat).subcortical_map_old])+10e-10];
    end

    %corticals
    g= figure,plot_cortical(HRF_map_cortical, 'surface_name', 'conte69','cmap','RdBu_r','color_range',range_cort)
    sgtitle([labels_hrf(idx_feat) ' old'],'Interpreter','none')

    %subcorticals
    g_sub=figure, plot_subcortical_AAL3(HRF_maps(idx_feat).subcortical_map_old,'ventricles','False', 'cmap', 'RdBu_r','color_range',range_subcort)
    sgtitle([labels_hrf(idx_feat) ' old'],'Interpreter','none')

end

