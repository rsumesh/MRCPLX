%% Script to benchmark T2 fit execution times from simulated 8-channel T2* data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

warning on

load('../data/t2s_mrx_benchmarking_data.mat')
% Contains simulated data from different SNR and parameter combinations.
% Sub-sample data. We found that the execution times differs with data
% size. This may depend on Matlab options, computer ram and other
% parameters beyond what we could take care of. We sub-sampled the data by
% a factor of 3.
tes=tes(1:3:end);
t2s_mrx_data=t2s_mrx_data(:,1:3:end,:);

algo_names={'LC','VC'};
fit_routines=1:2;

exec_times=zeros(size(t2s_mrx_data,1),length(fit_routines));
t2s_est=zeros(size(t2s_mrx_data,1),length(fit_routines));

disp(['Total time series to fit: ' num2str(size(t2s_mrx_data,1))]);

for idx=1:size(t2s_mrx_data,1)
    disp(num2str(idx))
    for routine_idx=1:length(fit_routines)
        [temp_exec_time,temp_t2s_est]=benchmark_t2s_mrx_exec_time(t2s_mrx_data(idx,:,:),tes,fit_routines(routine_idx));
        exec_times(idx,routine_idx)=temp_exec_time;
        t2s_est(idx,routine_idx)=temp_t2s_est;
    end
end

t2s_est_errors=t2s_est-repmat(T2s,[1 length(fit_routines)]);
t2s_est_errors=100*t2s_est_errors./repmat(T2s,[1 length(fit_routines)]);% In percentage
save('t2s_benchmarking_results_mrx.mat','exec_times','t2s_est','T2s','t2s_est_errors');