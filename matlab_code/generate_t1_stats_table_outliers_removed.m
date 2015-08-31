function outliers=generate_t1_stats_table_outliers_removed(all_algos,fp,r1_fit_results,snr,init_phase,T1,algo,weight_flag,noise_thresh_flag)
% generate_t1_stats_table_outliers_removed: Function to remove outlies and write the results to a csv file
% usage: outliers=generate_t1_stats_table_outliers_removed(all_algos,fp,r1_fit_results,snr,init_phase,T1,algo,weight_flag,noise_thresh_flag)
% Inputs:
% all_algos: List of all algorithms, fp: file pointer (to the csv file)
% r1_fit_results: Structure returned by calculate_r1_all_methods function
% snr: simulation SNR, init_phase: True initial phase (degrees), T1: True T1
% algo: A number identifying the algorithm used to fit the data
% weight_flag: Determines if weighted-fits were performed (1) or not (0)
% noise_thresh_flag: Determines if noise thresholding was done (1) or not (0)
% Output:
% outliers: A structure containing
% indices: Locations where the outliers exist in the r1_fit_results
% vals: Outlier values
% out_frac: Fraction of outliers (from all the fit results)

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

start_ph=round(init_phase*180/pi);

if ~isfield(r1_fit_results,'ph')%Magnitude fit
    vals=cat(2,1./r1_fit_results.r1,r1_fit_results.m0,r1_fit_results.bbya);
    params={'t1est','m0est','bbya'};
    true_params=[T1;1.0;2.0];% Proton density was set to 1 and parameter p (bbya) was set to 2 in simulations
else%Complex-fit
    vals=cat(2,1./r1_fit_results.r1,r1_fit_results.bbya,(180/pi)*r1_fit_results.ph,r1_fit_results.m0);
    params={'t1est','bbya','phest','m0est'};
    true_params=[T1;2.0;start_ph;1.0];
end

% Find which of the estimates are outliers
[~,temp_indices,~]=deleteoutliers_iqr(1./r1_fit_results.r1);
outliers.indices=temp_indices;% Location of outliers
outliers.vals=single(vals(temp_indices,:));% Values of outliers
outliers.out_frac=length(temp_indices)/size(vals,1);% Fraction of outliers

disp(['Outlier_fraction: ' num2str(outliers.out_frac)]);


% Statements to write parameter estimates from all algorithms to a single
% file. Comment them if you do not want this feature.
weight_str=num2str(weight_flag);
nt_str=num2str(noise_thresh_flag);
algo_name=all_algos{algo};
t1_str=num2str(T1);
ph_str=num2str(start_ph);
snr_str=num2str(snr);
vox_str=cell(1,size(r1_fit_results.r1,1));
for idx=1:length(vox_str)
    vox_str{idx}=num2str(idx);
end
vals(temp_indices,:)=[];

for vox_idx=1:size(vals,1)
    for param_idx=1:size(vals,2)
        fprintf(fp,'%s;',vox_str{vox_idx});
        fprintf(fp,'%s;',algo_name);
        fprintf(fp,'%s;',t1_str);
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