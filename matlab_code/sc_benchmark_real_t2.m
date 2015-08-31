%% Script to benchmark T2 fitting times of different algorithms
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

load('../data/t2_benchmarking_data.mat');%Load benchmarking data (1000 example time series from SNR-degraded real data)

all_algos={'LG','LQ','AR','NM','LC','NC','VC','LP'};% List of all fitting algorithms used, see Table 1 of the manuscript for details. We have omitted RM since our implementation was too slow.
fit_routines=1:8;%Pick the methods that you choose to fit the data with.

exec_times=zeros(size(t2_data,1),length(fit_routines));
t2_est=exec_times;

warning off
for idx=1:size(t2_data,1)
    disp(num2str(idx));
    for routine_idx=1:length(fit_routines)
        % Function that was used to obtain the benchmarking results
        [temp_exec_time,temp_t2_est]=benchmark_t2_exec_time(t2_data(idx,:),te,fit_routines(routine_idx));
        exec_times(idx,routine_idx)=temp_exec_time;
        t2_est(idx,routine_idx)=temp_t2_est;
    end
end

% Comment or modify the save statement suitably
save('t2_benchmarking_results.mat','exec_times','t2_est');