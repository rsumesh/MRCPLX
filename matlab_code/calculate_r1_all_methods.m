function est_results=calculate_r1_all_methods(data,ti,fit_type,noise_std,weight_flag,noise_thresh_flag,varargin)
%CALCULATE_R1_ALL_METHODS: Methods to do R1 fits (magnitude-based and complex-valued data-based)
%usage: est_results=CALCULATE_R1_ALL_METHODS(data,ti,fit_type,noise_std,weight_flag,noise_thresh_flag,varargin)
% Inputs:
% data: In the format (x,y,t) where t denotes the temporal dimension (different inversion time points)
% ti: inversion time points (s)
% fit_type: A number indicating the fitting method (1-10)
% All methods ending in MAG denote that they work with magnitude data
% All methods ending in CPLX denote that they work with complex-valued data
% Options: 1: LQ-MAG Levenberg-Marquardt fit, 2: NM-MAG NIPALS Gaussian noise assumption,
% 3: RM-MAG NIPALS Rician noise assumption, 4: LS1-CPLX Derivative-free Least-squares type 1 (not discussed in the article)
% 5: LC2-CPLX Derivative-based Least-squares type 2 (not discussed in the article)
% 6: LC-CPLX Derivative-based Least-squares, 7: NC-CPLX NIPALS on complex-valued data
% 8: VC-CPLX1 VARPRO Type 1 on complex-valued data (not discussed in the article)
% 9: VC-CPLX2 VARPRO Type 2 on complex-valued data (not discussed in the article)
% 10: VC-CPLX VARPRO on complex-valued data
% noise_std: standard deviation of the noise (common for both real and
% imaginary channels)
% weight_flag: Whether to do weighted estimates or not
% noise_thresh_flag: Whether to do noise-thresholding or not
% varargin: If provided, a mask (helpful for 3-d (x,y,t) datasets)
%
% Output:
% est_results: A structure containing the estimated parameters and
% timing information. Structure fields can be
% r1: estimated r1 values, m0: estimated proton density,
% rel_res: weighted residual norms (for VARPRO-based fits)
% eflag: exit flags (when optimizers are used)
% fval: objective value (esimated relative residue (2-norm)) at the end of
% optimization (if optimizers are involved)
% npts_fit: No. of points used for the fit (useful in case noise-thresholding was
% invoked)
% ph: Estimated initial phase (degrees)
% tTimeit: Estimated execution time using timeit command on the first time series.
% tictoc: Estimated execution time using tic and toc commands (for fitting all the time series)
% tcpu: Estimated cpu execution time using cputime command (for fitting all the time series)

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

r1_lb=0.1;
r1_ub=10;
ti=ti(:);
noise_thresh=3*noise_std;

if ~isempty(varargin)
    mask=varargin{1,1};
else
    mask=ones(size(data,1),size(data,2));% If no mask is provided, fit all time series
end

% Optional statements to run the profiler
% profile off;
% profile clear;
% profile on;

tStart = tic;% Start clock to estimate time to fit all the time series
tStartcpu = cputime;
tic

switch fit_type
    case 1
        disp('LQ-MAG Levenberg-Marquardt fit');
        [r1_fit_results,tTimeit]=calculate_r1_LQ(data,ti,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 2
        disp('NM-MAG NIPALS Gaussian Assumption');
        [r1_fit_results,tTimeit]=calculate_r1_NM(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 3
        disp('RM-MAG NIPALS Rician assumption');
        [r1_fit_results,tTimeit]=calculate_r1_RM(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,noise_std,mask);
        
    case 4
        disp('LC1-CPLX Derivative-free Least-squares');%Not discussed in the article
        [r1_fit_results,tTimeit]=calculate_r1_LC1(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 5
        disp('LC2-CPLX Derivative-based Least-squares');%Not discussed in the article
        [r1_fit_results,tTimeit]=calculate_r1_LC2(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 6
        disp('LC-CPLX Derivative-based Least-squares');
        [r1_fit_results,tTimeit]=calculate_r1_LC(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 7
        disp('NC-CPLX NIPALS on complex-valued data');
        [r1_fit_results,tTimeit]=calculate_r1_NC(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 8
        disp('VC-CPLX1 VARPRO Type 1 on complex-valued data');%Not discussed in the article
        [r1_fit_results,tTimeit]=calculate_r1_VC1(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 9
        disp('VC-CPLX2 VARPRO Type 2 on complex-valued data');%Not discussed in the article
        [r1_fit_results,tTimeit]=calculate_r1_VC2(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    case 10
        disp('VC-CPLX VARPRO on complex-valued data');
        [r1_fit_results,tTimeit]=calculate_r1_VC(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,mask);
        
    otherwise
        disp('Unknown Algorithm!');
end

r1_fit_results(tTimeit==Inf,:,:,:,:)=[];

if size(r1_fit_results,3)==7 % Fits to magnitude data
    est_results.r1=r1_fit_results(:,:,1);
    est_results.bbya=r1_fit_results(:,:,2);
    est_results.m0=r1_fit_results(:,:,3);
    est_results.rel_res=r1_fit_results(:,:,4);
    est_results.eflag=r1_fit_results(:,:,5);
    est_results.fval=r1_fit_results(:,:,6);
    est_results.npts_fit=r1_fit_results(:,:,7);
elseif size(r1_fit_results,3)==8 % Fits to complex-valued data
    est_results.r1=r1_fit_results(:,:,1);
    est_results.bbya=r1_fit_results(:,:,2);
    est_results.ph=r1_fit_results(:,:,3);
    est_results.m0=r1_fit_results(:,:,4);
    est_results.rel_res=r1_fit_results(:,:,5);
    est_results.eflag=r1_fit_results(:,:,6);
    est_results.fval=r1_fit_results(:,:,7);
    est_results.npts_fit=r1_fit_results(:,:,8);
end

toc
tElapsed = toc(tStart);% Calculate the time taken for the fitting to complete
tElapsedcpu=cputime -tStartcpu;

est_results.tTimeit=tTimeit;
est_results.tictoc=tElapsed;
est_results.tcpu=tElapsedcpu;

% Optional statements to view the profiler output
% profile off
% profile viewer

temp=est_results.r1;
temp=temp(:);

if ~isempty(mask)
    temp=temp(logical(mask(:)));
end

% Display some results
disp(['Median R1 (1/s): ' num2str(median(temp))]);
disp(['Median T1 (ms): ' num2str(1000/median(temp))]);

disp(['tictoc time (s): ' num2str(est_results.tictoc)]);
disp(['cpu time (s): ' num2str(est_results.tcpu)]);
disp(['timeit time (ms): ' num2str(1000*tTimeit)]);

return;