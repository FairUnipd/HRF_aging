function [params_y,params_o,HRFmatrix_y,HRFmatrix_o,Y_LR_y,Y_LR_o] = HRF_partition(files,age_young,age_old)
%performs HRF partition into young and old subsets

    county=1;
    counto=1;
    countall=1;
    
    for sub=1:length(files)
      
       load(files(sub).name,'age_sub','hrf_params','h')
    
       %exclude middle-aged subjects
       if age_sub>age_young && age_sub<age_old
           continue;
       end
    
       for ii=1:size(hrf_params,2)
    
           for jj=1:length(node_bounds)
                if age_sub<=age_young
                    params_y(ii,jj,county)=(hrf_params(node_bounds(jj),ii));
                    HRFmatrix_y(:,jj,county)=(h(:,node_bounds(jj)));
                    
                end
                if age_sub>=age_old
                    params_o(ii,jj,counto)=(hrf_params(node_bounds(jj),ii));
                    HRFmatrix_o(:,jj,counto)=(h(:,node_bounds(jj)));
                end
           end
           
       end
       %partition age value
       if age_sub<=age_young
                Y_LR_y(county)=age_sub;
                county=county+1;
       end
       if age_sub>=age_old
                Y_LR_o(counto)=age_sub;  
                counto=counto+1;
       end
    
       clear age_sub hrf_params
    
    end
end