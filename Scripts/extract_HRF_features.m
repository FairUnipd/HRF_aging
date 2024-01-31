%% HRF params
%example script to extract subject-level the 9 HRF features 
%h = matrix of HRF curves from sDCM (num of time points x num of ROIs)
%TR = repetition time

addpath('Data')
load example_subj.mat
%Peak amplitude
[hrf_Ap,hrf_idx]=max(h);

%time to peak
hrf_tp=hrf_idx*TR;

for ii=1:num_rois %for each ROI

    %trough amplitude
    [maxA,idx_max]=max(h(:,ii));
    minA=min(h(idx_max:end-1,ii));
    hrf_At(ii)=minA;
    idx_min=find(h(:,ii)==minA);
    
    %time to trough
    hrf_tt(ii)=idx_min*TR;

    %peak to trough
    hrf_p2t(ii)=hrf_tt(ii)-hrf_tp(ii);

    %rate features: approximate HRF curve by a triangle 
    
    %fall slope (slope between max value and value computed after 3 TR,
    %returns the best results by visual inspection
    p_fall=polyfit([idx_max*TR,(idx_max+3)*TR],[maxA,h(idx_max+3,ii)],1);
    hrf_fall_slope(ii)=p_fall(1);

    %rise slope
    hrf_diff=diff(h(1:idx_max,ii));
    idx_first_slope=find(hrf_diff<0,1,'last');
    if ~isempty(idx_first_slope)
        x_ax=[(idx_first_slope+1)*TR,idx_max*TR];%(idx_first_slope+1)*TR:TR:idx_max*TR;
        p_rise=polyfit([(idx_first_slope+1)*TR,idx_max*TR],[h(idx_first_slope+1,ii),maxA],1);
    else
        x_ax=[TR,idx_max*TR];%TR:TR:idx_max*TR;
        p_rise=polyfit([TR,idx_max*TR],[h(1,ii),maxA],1);  
    end
    hrf_rise_slope(ii)=p_rise(1);

    %AUC
    hrf_auc(ii)=trapz(TR,h(:,ii));

    %FWHM
    [pks,locs,w,p] = findpeaks(h(:,ii),1/TR);
    idx_first_peak=find(pks==maxA);
    hrf_fwhm(ii)=w(idx_first_peak);

%rate features sanity check
%     figure
%     plot(TR:TR:TR*length(h(:,ii)),h(:,ii))
%     hold on
%     plot(x_ax,polyval(p_rise,x_ax))
%     hold on
%     plot([idx_max*TR,idx_min*TR],polyval(p_fall,[idx_max*TR,idx_min*TR]))
%     hold off

%AUC sanity check
%     figure
%     area(TR:TR:TR*length(h(:,ii)),h(:,ii))
end

%save HRF features (at subject level)
hrf_params=[hrf_Ap;hrf_At;hrf_tp;hrf_tt;hrf_rise_slope;hrf_fall_slope;...
hrf_auc;hrf_fwhm;hrf_p2t]';