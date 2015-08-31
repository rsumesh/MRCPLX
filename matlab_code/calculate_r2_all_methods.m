function est_results=calculate_r2_all_methods(data,tau,fit_type,noise_std,max_f_deviation,weight_flag,noise_thresh_flag,varargin)
%CALCULATE_R2_ALL_METHODS: Methods to do R2* (and R2) fits (magnitude-based
%and complex-valued data-based)
%usage: est_results=CALCULATE_R2_ALL_METHODS(data,tau,fit_type,noise_std,max_f_deviation,weight_flag,noise_thresh_flag,varargin)
% Inputs:
% data: In the format (x,y,t) where t denotes the temporal dimension (different echoes)
% tau: echo_spacing (s)
% fit_type: A number indicating the fitting method (1-11)
% All methods ending in MAG denote that they work with magnitude data
% All methods ending in CPLX denote that they work with complex-valued data
% Options: 1: LG-MAG Linear fit using logarithm, 2: LQ-MAG Levenberg-Marquardt fit
% 3: AR-MAG ARLO fit, 4: NM-MAG NIPALS Gaussian noise assumption,
% 5: VM-MAG VARPRO Gaussian assumption, 6: RM-MAG NIPALS Rician noise assumption
% 7: LS1-CPLX Derivative-free Least squares (not discussed in the article)
% 8: LC-CPLX Derivative-based Least squares
% 9: NC-CPLX NIPALS fit to complex-valued data (not discussed in the article)
% 10: VC-CPLX VARPRO fit to complex-valued data
% 11: LP-CPLX Linear prediction using Steiglitz-McBride iteration
% See Reference (below) for further details
% noise_std: standard deviation of the noise (common for both real and
% imaginary channels)
% max_f_deviation: Maximum frequency deviation that can be estimated (from
% complex-valued data)
% weight_flag: Whether to do weighted estimates or not
% noise_thresh_flag: Whether to do noise-thresholding or not
% varargin: If provided, a mask (helpful for 3-d (x,y,t) datasets)
%
%
% Output:
% est_results: A structure containing the estimated parameters and
% timing information. Structure fields can be
% r2: estimated r2 values, m0: estimated proton density,
% rel_res: weighted residual norms (for VARPRO-based fits)
% eflag: exit flags (when optimizers are used)
% fval: objective value (esimated relative residue (2-norm)) at the end of
% optimization (if optimizers are involved)
% npts_fit: No. of points used for the fit (useful in case noise-thresholding was
% invoked)
% ph: Estimated initial phase (degrees)
% fdev: Estimated frequency (B0) deviation (Hz)
% tTimeit: Estimated execution time using timeit command on the first time series.
% tictoc: Estimated execution time using tic and toc commands (for fitting all the time series)
% tcpu: Estimated cpu execution time using cputime command (for fitting all the time series)

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

noise_thresh=3*noise_std;% Threshold to discard data in case of noise-thresholding
r2_upper_threshold=200;% We do not expect T2* (T2) values below 5 ms (not in the simulations at least)

szdata=size(data);
te=tau*(1:szdata(3));

if isempty(max_f_deviation)
    max_f_deviation=0.5/max(diff(te(:)));
end

% Optional statements to run the profiler
% profile off
% profile clear
% profile on

tStart = tic; % Start clock to estimate time to fit all the time series
tStartcpu = cputime;

if ~isempty(varargin)
    mask=varargin{1};
else
    mask=ones(size(data,1),size(data,2));% If no mask is provided, fit all time series
end

switch fit_type
    case 1
        disp('LG-MAG Linear fit using logarithm');
        [r2_fit_results,tTimeit]=calculate_r2_LG(data,te,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 2
        disp('LQ-MAG Levenberg-Marquardt fit');
        [r2_fit_results,tTimeit]=calculate_r2_LQ(data,te,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 3
        disp('AR-MAG ARLO fit');
        [r2_fit_results,tTimeit]=calculate_r2_AR(data,te,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 4
        disp('NM-MAG NIPALS Gaussian Assumption');
        [r2_fit_results,tTimeit]=calculate_r2_NM(data,te,r2_upper_threshold,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 5
        disp('VM-MAG VARPRO Gaussian assumption');%Not discussed in the article
        [r2_fit_results,tTimeit]=calculate_r2_VM(data,te,r2_upper_threshold,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 6
        disp('RM-MAG NIPALS Rician assumption');
        [r2_fit_results,tTimeit]=calculate_r2_RM(data,te,r2_upper_threshold,weight_flag,noise_thresh_flag,noise_thresh,noise_std,mask);
        r2_fit_results=single(r2_fit_results);
        
    case 7
        disp('LS1-CPLX Derivative-free Least squares');%Not discussed in the article
        [r2_fit_results,tTimeit]=calculate_r2_LC1(data,te,r2_upper_threshold,max_f_deviation,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 8
        disp('LC-CPLX Derivative-based Least squares');
        [r2_fit_results,tTimeit]=calculate_r2_LC(data,te,r2_upper_threshold,max_f_deviation,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 9
        disp('NC-CPLX NIPALS fit to complex-valued data');
        [r2_fit_results,tTimeit]=calculate_r2_NC(data,te,r2_upper_threshold,max_f_deviation,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 10
        disp('VC-CPLX VARPRO fit to complex-valued data');
        [r2_fit_results,tTimeit]=calculate_r2_VC(data,te,r2_upper_threshold,max_f_deviation,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 11
        disp('LP-CPLX Linear prediction using Steiglitz-McBride iteration');
        [r2_fit_results,tTimeit]=calculate_r2_LP(data,te,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    otherwise
        disp('Unknown Algorithm!');
end


if size(r2_fit_results,3)==6 % Magnitude fits
    est_results.r2=r2_fit_results(:,:,1);
    est_results.m0=r2_fit_results(:,:,2);
    est_results.rel_res=r2_fit_results(:,:,3);
    est_results.eflag=r2_fit_results(:,:,4);
    est_results.fval=r2_fit_results(:,:,5);
    est_results.npts_fit=r2_fit_results(:,:,6);
    
elseif size(r2_fit_results,3)==8 % Complex-valued data fits
    est_results.r2=r2_fit_results(:,:,1);
    est_results.ph=r2_fit_results(:,:,2);
    est_results.fdev=r2_fit_results(:,:,3);
    est_results.m0=r2_fit_results(:,:,4);
    est_results.rel_res=r2_fit_results(:,:,5);
    est_results.eflag=r2_fit_results(:,:,6);
    est_results.fval=r2_fit_results(:,:,7);
    est_results.npts_fit=r2_fit_results(:,:,8);
else
    disp('Unknown output format!')
    return;
end

tElapsed = toc(tStart);% Calculate the time taken for the fitting to complete
tElapsedcpu=cputime -tStartcpu;

est_results.tTimeit=tTimeit;
est_results.tictoc=tElapsed;
est_results.tcpu=tElapsedcpu;

% Optional statements to view the profiler output
% profile off
% profile viewer

temp_result=est_results.r2;% For display purposes, we are only interested in R2* (R2) or equivalently T2* (T2)
temp_result=temp_result(:);
temp_result=temp_result(logical(mask(:)));

% Display some results
disp(['Median R2 (1/s): ' num2str(median(temp_result))]);
disp(['Median T2 (ms): ' num2str(1000/median(temp_result))]);

disp(['tictoc time (s): ' num2str(est_results.tictoc)]);
disp(['cpu time (s): ' num2str(est_results.tcpu)]);
disp(['timeit time (ms): ' num2str(1000*tTimeit)]);

return;