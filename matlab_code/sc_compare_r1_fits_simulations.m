%% Script to perform T1 fitting using various (magnitude and complex-valued data-based) methods

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

t1_data_path='../data/t1_data_sim.mat';%Path to the mat-file containing T1 data
op_folder='';%Directory where the outputs will be saved
op_filename='stats_table_r1_fit_results';%A csv file that contains parameter estimates from different fitting routines. Can be read in R, xl or any other spread-sheet software.

load(t1_data_path);

all_algos={'LQ','NM','RM','LC1','LC2','LC','NC','VC1','VC2','VC'};% List of all fitting algorithms used, see Table 1 of the manuscript for details.
% LC2, VC1 and VC2 are different variants of LC and VC, respectively, and have not
% been discussed in the manuscript. Their performance was the same as LC
% and VC. LC1 has also been left out from the manuscript results.

fit_routines=[1:3 6 7 10];%Pick the methods that you choose to fit the data with. We have omitted LC1, LC2, VC1 and VC2 from the list.

weight_flag=0;% For weighted fitting: 0-No weighting, 1-with weighting. Can also be arrayed, e.g. [0 1].
noise_thresh_flag=0;% For fits with noise-thresholding: 0-No noise thresholding, 1-with noise thresholding (3*sigma is the cut-off magnitude). Can be arrayed, e.g. [0 1].

% Create op_folder if it does not exist
if ~isempty(op_folder)
    if ~exist(op_folder,'dir')
        command_str=['mkdir -p ' op_folder];
        [status,result]=unix(command_str);
        if status
            disp(result);
        end
    end
end

fp=fopen([op_folder op_filename '.csv'],'w');%File that contains all the fitting results

header={'vox','algo','r1','ph','wf','ntf','snr','param','evalue','tvalue'};% The header of the csv file
% vox: voxel number, algo: Algorithm, r1: True R1 value,
% ph: starting phase, wf: weighting flag, ntf: noise threshold flag, snr: signal to noise ratio,
% param: One of estimated R1, proton density or starting phase and more
% evalue: estimated value of the parameter, tvale: true value of the parameter

% Write the header (first line) to the csv file
for fields_idx=1:length(header)
    if fields_idx~=length(header)
        fprintf(fp,'%s;',header{fields_idx});
    else
        fprintf(fp,'%s\n',header{fields_idx});
    end
end

szdata=size(t1_data);

warning on
useful_fields={'r1','m0','bbya','ph','tTimeit','tictoc','tcpu'};% Parameters of interest

r1_fit_results=cell([length(fit_routines) szdata([2 4:5]) length(weight_flag) length(noise_thresh_flag)]);% Create a cell to hold all the results (in case we need to save it to a matfile)
outliers=r1_fit_results;% Create a cell to hold all outliers as well

for algo_idx=1:length(fit_routines)% For each selected algorithm
    for snr_idx=1:szdata(5)% For each selected SNR
        for ph_idx=1:szdata(4)% For every starting phase condition
            for r1_idx=1:szdata(2)% For every simulated T1 (R1)
                for noise_thresh_idx=1:length(noise_thresh_flag)% For different noise-thresholding conditions
                    for weight_idx=1:length(weight_flag)% For different weighting conditions
                        
                        disp_str1=['SNR: ' num2str(snrs(snr_idx))];
                        disp_str3=['ph: ' num2str(start_ph(ph_idx))];
                        disp_str4=['T1 (ms): ' num2str(1e3*T1(r1_idx))];
                        disp_str5=['R1 (1/s): ' num2str(1/T1(r1_idx))];
                        disp_str6=['Algo: ' all_algos{fit_routines(algo_idx)}];
                        disp_str7=['Noise-thresh-flag: ' num2str(noise_thresh_flag(noise_thresh_idx))];
                        disp_str8=['Weight-flag: ' num2str(weight_flag(weight_idx))];
                        
                        disp_str=[disp_str6 ' ' disp_str1 ' ' disp_str4 ' ' disp_str5 ' ' disp_str3 ' ' disp_str7 ' ' disp_str8];
                        disp(disp_str);% Display the fitting condition
                        temp_data=squeeze(t1_data(:,r1_idx,:,ph_idx,snr_idx));
                        temp_data=reshape(temp_data,[szdata(1) 1 szdata(3)]);
                        
                        % The main function which contains all other fitting-related information
                        temp_results=calculate_r1_all_methods...
                            (temp_data,ti(:),fit_routines(algo_idx),noise_std(snr_idx),weight_flag(weight_idx),noise_thresh_flag(noise_thresh_idx));
                        
                        temp_field_names=fieldnames(temp_results);
                        unwanted_fields = temp_field_names(~ismember(temp_field_names,useful_fields));
                        temp_results=rmfield(temp_results,unwanted_fields);
                        
                        % Function that weeds out the outliers and also
                        % writes the useful results to the csv file
                        temp_outliers=generate_t1_stats_table_outliers_removed(all_algos,fp,temp_results,snrs(snr_idx),start_ph(ph_idx),T1(r1_idx),fit_routines(algo_idx),weight_flag(weight_idx),noise_thresh_flag(noise_thresh_idx));
                        disp(' ');%For pretty printing
                        
                        r1_fit_results{algo_idx,r1_idx,ph_idx,snr_idx,weight_idx,noise_thresh_idx}=temp_results;
                        outliers{algo_idx,r1_idx,ph_idx,snr_idx,weight_idx,noise_thresh_idx}=temp_outliers;
                    end
                end
            end
        end
    end
end
fclose(fp);

% Optional statements to save the results to mat files (beware, these can grow big very quickly)
save([op_folder 'r1_fit_results.mat'],'r1_fit_results');
save([op_folder 'r1_outlier_results.mat'],'outliers');
