%% Script to perform T2* fitting using various (magnitude and complex-valued data-based) methods

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

t2s_data_path='../data/t2s_data_sim.mat';%Path to the mat-file containing t2s data
op_folder='';%Directory where the outputs will be stored
op_filename='stats_table_r2s_fit_results';%A csv file that contains parameter estimates from different fitting routines. Can be read in R, xl or any other spread-sheet software.

load(t2s_data_path);

all_algos={'LG','LQ','AR','NM','VM','RM','LS1','LC','NC','VC','LP'};% List of all fitting algorithms used, see Table 1 of the manuscript for details.
% VM and LS1 are not discussed in the manuscript as they seemed redundant.

fit_routines=[1:4 6 8:11];%Pick the methods that you choose to fit the data with. We have omitted VM and LS1 from the list.

weight_flag=0;% For weighted fitting: 0-No weighting, 1-with weighting. Can also be arrayed, e.g. [0 1].
noise_thresh_flag=0;% For fits with noise-thresholding: 0-No noise thresholding, 1-with noise thresholding (3*sigma is the cut-off magnitude). Can be arrayed, e.g. [0 1].
max_freq_deviation=0.5/te;%Hz, maximum frequency deviation that can be estimated with complex-valued data (Nyquist frequency)

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

header={'vox','algo','r2s','foffset','ph','wf','ntf','snr','param','evalue','tvalue'};% The header of the csv file
% vox: voxel number, algo: Algorithm, r2s: True R2* value,
% foffset: Frequency deviation, ph: starting phase,
% wf: weighting flag, ntf: noise threshold flag, snr: signal to noise ratio,
% param: One of estimated R2*, proton density or frequency deviation or starting phase and more
% evalue: estimated value of the parameter, tvale: true value of the parameter

% Write the header (first line) to the csv file
for fields_idx=1:length(header)
    if fields_idx~=length(header)
        fprintf(fp,'%s;',header{fields_idx});
    else
        fprintf(fp,'%s\n',header{fields_idx});
    end
end

szdata=size(t2s_data);
R2s=1./T2s;

warning on
useful_fields={'r2','m0','ph','fdev','tTimeit','tictoc','tcpu'};% Parameters that would be written to the the csv file if they exist

r2s_fit_results=cell([length(fit_routines) szdata([2 4:6]) length(weight_flag) length(noise_thresh_flag)]);% Create a cell to hold all the results (in case we need to save it to a matfile)
outliers=r2s_fit_results;% Create a cell to hold all outliers as well

for algo_idx=1:length(fit_routines)% For each selected algorithm
    for snr_idx=1:szdata(6)% For each selected SNR
        for f_idx=1:szdata(4)% For each frequency deviation condition
            for ph_idx=1:szdata(5)% For every starting phase condition
                for r2s_idx=1:szdata(2)% For every simulated T2* (R2*)
                    for noise_thresh_idx=1:length(noise_thresh_flag) % For different noise-thresholding conditions
                        for weight_idx=1:length(weight_flag)% For different weighting conditions
                            
                            disp_str1=['SNR: ' num2str(snrs(snr_idx))];
                            disp_str2=['Freq-offset(Hz): ' num2str(f_range(f_idx))];
                            disp_str3=['ph: ' num2str(start_ph(ph_idx))];
                            disp_str4=['T2s (ms): ' num2str(1e3*T2s(r2s_idx))];
                            disp_str5=['R2s (1/s): ' num2str(1/T2s(r2s_idx))];
                            disp_str6=['Algo: ' all_algos{fit_routines(algo_idx)}];
                            disp_str7=['Noise-thresh-flag: ' num2str(noise_thresh_flag(noise_thresh_idx))];
                            disp_str8=['Weight-flag: ' num2str(weight_flag(weight_idx))];

                            disp_str=[disp_str6 ' ' disp_str1 ' ' disp_str4 ' ' disp_str5 ' ' disp_str2 ' ' disp_str3 ' ' disp_str7 ' ' disp_str8];
                            disp(disp_str); % Display the fitting condition
                            
                            temp_data=squeeze(t2s_data(:,r2s_idx,:,f_idx,ph_idx,snr_idx)); % Select the appropriate data
                            temp_data=reshape(temp_data,[szdata(1) 1 szdata(3)]);% Reshape it to the format the fitting algorithms will understand it
                            
                            % The main function which contains all other fitting-related information
                            temp_results=calculate_r2_all_methods...
                                (temp_data,te,fit_routines(algo_idx),noise_std(snr_idx),max_freq_deviation,weight_flag(weight_idx),noise_thresh_flag(noise_thresh_idx));
                            
                            temp_field_names=fieldnames(temp_results);% See what results are returned, not all routines can return all the fields.
                            unwanted_fields = temp_field_names(~ismember(temp_field_names,useful_fields));
                            temp_results=rmfield(temp_results,unwanted_fields);% Retain only useful fields
                            
                            % Function that weeds out the outliers and also
                            % writes the useful results to the csv file
                            temp_outliers=generate_t2s_stats_table_outliers_removed(all_algos,fp,temp_results,snrs(snr_idx),f_range(f_idx),start_ph(ph_idx),T2s(r2s_idx),fit_routines(algo_idx),weight_flag(weight_idx),noise_thresh_flag(noise_thresh_idx));
                            disp(' ');
                            
                            r2s_fit_results{algo_idx,r2s_idx,f_idx,ph_idx,snr_idx,weight_idx,noise_thresh_idx}=temp_results;% Save the results in a cell
                            outliers{algo_idx,r2s_idx,f_idx,ph_idx,snr_idx,weight_idx,noise_thresh_idx}=temp_outliers;% Save the outliers in a cell
                        end
                    end
                end
            end
        end
    end
end
fclose(fp);

% Optional statements to save the results to mat files (beware, these can grow big very quickly)
save([op_folder 'r2s_fit_results.mat'],'r2s_fit_results');
save([op_folder 'r2s_outlier_results.mat'],'outliers');
