%% Script to benchmark T1 fitting times of different algorithms
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

load('../data/t1_benchmarking_data.mat');%Load benchmarking data (1000 example time series from SNR-degraded real data)

all_algos={'LQ','NM','LC','NC','VC'};% List of all fitting algorithms used, see Table 1 of the manuscript for details. We have omitted RM since our implementation was too slow.
fit_routines=1:5;%Pick the methods that you choose to fit the data with.

exec_times=zeros(size(t1_data,1),length(fit_routines));
t1_est=exec_times;

for idx=1:size(t1_data,1)
    disp(num2str(idx));
    for routine_idx=1:length(fit_routines)
        % Function that was used to obtain the benchmarking results
        [temp_exec_time,temp_t1_est]=benchmark_t1_exec_time(t1_data(idx,:),ti,fit_routines(routine_idx));
        exec_times(idx,routine_idx)=temp_exec_time;
        t1_est(idx,routine_idx)=temp_t1_est;
    end
end

% Comment or modify the save statement suitably
save('t1_benchmarking_results.mat','exec_times','t1_est');