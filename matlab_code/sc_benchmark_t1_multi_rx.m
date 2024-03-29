%% Script to benchmark T1 fit execution times from simulated 8-channel T1 data

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

warning on
load('../data/t1_mrx_benchmarking_data.mat');

fit_routines=1:2;
algo_names={'LC','VC'};

% We found that the execution times differs with data
% size. This may depend on Matlab options, computer ram and other
% parameters beyond what we could account for. Please be aware of this
% while comparing the results to those reported in the article.
mrx_t1_data(:,2:2:end,:)=[];
mrx_t1_data(:,2:2:end,:)=[];
mrx_t1_data(:,2:2:end,:)=[];

ti(2:2:end)=[];
ti(2:2:end)=[];
ti(2:2:end)=[];

exec_times=zeros(size(mrx_t1_data,1),length(fit_routines));
t1_est=zeros(size(mrx_t1_data,1),length(fit_routines));

disp(['Num. time series to fit: ' num2str(size(mrx_t1_data,1))]);

for idx=1:size(mrx_t1_data,1)
    disp(num2str(idx))
    for routine_idx=1:length(fit_routines)
        [temp_exec_time,temp_t1_est]=benchmark_mrx_t1_exec_time(mrx_t1_data(idx,:,:),ti,fit_routines(routine_idx));
        exec_times(idx,routine_idx)=temp_exec_time;
        t1_est(idx,routine_idx)=temp_t1_est;
    end
end

t1_est_errors=t1_est-repmat(T1s,[1 length(fit_routines)]);
t1_est_errors=100*t1_est_errors./repmat(T1s,[1 length(fit_routines)]);% In percentage
save('t1_mrx_benchmarking_results.mat','exec_times','t1_est','T1s','t1_est_errors','algo_names');
