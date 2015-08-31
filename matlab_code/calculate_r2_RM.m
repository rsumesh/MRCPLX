function [r2_fit_results,tTimeit,fit_errors]=calculate_r2_RML(data,te,t2_upper_threshold,weight_flag,noise_thresh_flag,noise_thresh,noise_std,varargin)

% Function to calculate R2 maps using NIPALS on magnitude data but using
% Rician noise assumption
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

noise_var=noise_std^2;
NLPstart=100e-3;
INLB=0;
INUB=t2_upper_threshold;
tTimeit=0;

szdata=size(data);
r2_fit_results=zeros([szdata(1:2) 6]);
fit_errors=zeros(szdata);

if ~isempty(varargin{1})
    temp=varargin{1};
else
    temp=ones(szdata(1:2));
end

[rdout_pos,phenc_pos]=find(temp);
clear temp

if weight_flag
    weights=abs(data);
else
    w=ones(size(data,3),1);
end

options=optimset('fminsearch');
options=optimset(options,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off');

for idx=1:length(phenc_pos)
    time_series=abs(squeeze(data(rdout_pos(idx),phenc_pos(idx),:)));
    
    if weight_flag
        w=squeeze(weights(rdout_pos(idx),phenc_pos(idx),:));
    end
    
    if noise_thresh_flag
        w(time_series<noise_thresh)=0;
    end
    
    fit_numvals=sum(double(logical(w)));
    
    if fit_numvals>2
        
        if idx==1 % One type of benchmarking
            f=@()fminspleas_rician({@(param,xdata)exp(-xdata(:)*param)},NLPstart,te(:),time_series(:),INLB,INUB,w(:),noise_var,options);
            tTimeit=timeit(f);
        end
        
        [INLP,ILP,eflag,fval]=fminspleas_rician({@(param,xdata)exp(-xdata(:)*param)},NLPstart,te(:),time_series(:),INLB,INUB,w(:),noise_var,options);
        predicted_time_series=ILP*exp(-te(:)*INLP);
        
        fit_errors(rdout_pos(idx),phenc_pos(idx),:)=time_series(:)-predicted_time_series(:);
        error=norm(squeeze(fit_errors(rdout_pos(idx),phenc_pos(idx),:)))/norm(time_series);
        
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[INLP;ILP;error;eflag;fval;fit_numvals];
    else
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
    
end

r2_fit_results=single(r2_fit_results);

return;