function outliers=generate_t2s_stats_table_outliers_removed(all_algos,fp,r2s_fit_results,snr,foffset,init_phase,T2s,algo,weight_flag,noise_thresh_flag)
% generate_t2s_stats_table_outliers_removed: Function to remove outlies and write the results to a csv file
% usage: outliers=generate_t2s_stats_table_outliers_removed(all_algos,fp,r2s_fit_results,snr,foffset,init_phase,T2s,algo,weight_flag,noise_thresh_flag)
% Inputs:
% all_algos: List of all algorithms, fp: file pointer (to the csv file)
% r2s_fit_results: Structure returned by calculate_r2_all_methods function
% snr: simulation SNR, foffset: True frequency (B0) offset (Hz)
% init_phase: True initial phase (degrees), T2s: True T2* (T2)
% algo: A number identifying the algorithm used to fit the data
% weight_flag: Determines if weighted-fits were performed (1) or not (0)
% noise_thresh_flag: Determines if noise thresholding was done (1) or not (0)
% Output:
% outliers: A structure containing
% indices: Locations where the outliers exist in the r2s_fit_results
% vals: Outlier values
% out_frac: Fraction of outliers (from all the fit results)

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY


start_ph=round(init_phase*180/pi);

if ~isfield(r2s_fit_results,'ph')%Magnitude fit
    vals=cat(2,1./r2s_fit_results.r2,r2s_fit_results.m0);
    params={'t2sest','m0est'};
    true_params=[T2s;1.0];% Proton density was set to 1 in simulations
else %Complex-fit
    vals=cat(2,1./r2s_fit_results.r2,r2s_fit_results.fdev,(180/pi)*r2s_fit_results.ph,r2s_fit_results.m0);
    params={'t2sest','fdev','phest','m0est'};
    true_params=[T2s;foffset;start_ph;1.0];
end

% Find which of the estimates are outliers
[~,temp_indices]=deleteoutliers_iqr(1./r2s_fit_results.r2);
outliers.indices=temp_indices;% Location of outliers
outliers.vals=single(vals(temp_indices,:));% Values of outliers
outliers.out_frac=length(temp_indices)/size(vals,1);% Fraction of outliers

disp(['Outlier_fraction: ' num2str(outliers.out_frac)]);


% Statements to write parameter estimates from all algorithms to a single
% file. Comment them if you do not want this feature.
weight_str=num2str(weight_flag);
nt_str=num2str(noise_thresh_flag);
algo_name=all_algos{algo};
t2s_str=num2str(T2s);
f_offset_str=num2str(foffset);
ph_str=num2str(start_ph);
snr_str=num2str(snr);
vox_str=cell(1,size(r2s_fit_results.r2,1));
for idx=1:length(vox_str)
    vox_str{idx}=num2str(idx);
end

vals(temp_indices,:)=[];
for vox_idx=1:size(vals,1)
    for param_idx=1:size(vals,2)
        fprintf(fp,'%s;',vox_str{vox_idx});
        fprintf(fp,'%s;',algo_name);
        fprintf(fp,'%s;',t2s_str);
        fprintf(fp,'%s;',f_offset_str);
        fprintf(fp,'%s;',ph_str);
        fprintf(fp,'%s;',weight_str);
        fprintf(fp,'%s;',nt_str);
        fprintf(fp,'%s;',snr_str);
        fprintf(fp,'%s;',params{param_idx});
        fprintf(fp,'%f;',vals(vox_idx,param_idx));
        fprintf(fp,'%f\n',true_params(param_idx));
    end
end

return;