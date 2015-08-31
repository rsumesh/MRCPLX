function [r2_fit_results,tTimeit]=calculate_r2_LG(data,te,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R2 maps using logarithm of magnitude data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

tTimeit=0;
szdata=size(data);

if ~isempty(varargin{1})
    temp=varargin{1};
else
    temp=ones(szdata(1:2));
end

[rdout_pos,phenc_pos]=find(temp);
clear temp

H=[ones(szdata(3),1) -te(:)];
w=ones(size(data,3),1);
pinvH=pinv(H);

r2_fit_results=zeros([szdata(1:2) 6]);

for idx=1:length(phenc_pos)
    time_series=abs(squeeze(data(rdout_pos(idx),phenc_pos(idx),:)));
    log_time_series=log(time_series);
    
    if weight_flag
        w=time_series(:);
        if noise_thresh_flag
            w(w<noise_thresh)=0;
        end
        W_mtx=diag(w.^0.5);
        pinvH=pinv(W_mtx*H);
        log_time_series=W_mtx*log_time_series(:);
    elseif noise_thresh_flag
        w(time_series(:)<noise_thresh)=0;
        W_mtx=diag(w);
        pinvH=pinv(W_mtx*H);
        log_time_series=W_mtx*log_time_series(:);
    end
    
    if idx==1
        f=@()(pinvH*log_time_series(:));
        tTimeit=timeit(f);
    end
    
    optimal_params=pinvH*log_time_series(:);
    
    predicted_time_series=exp(optimal_params(1))*exp(-te(:)*optimal_params(2));
    
    r2_fit_results(rdout_pos(idx),phenc_pos(idx),1)=optimal_params(2);
    r2_fit_results(rdout_pos(idx),phenc_pos(idx),2)=exp(optimal_params(1));
    r2_fit_results(rdout_pos(idx),phenc_pos(idx),3)=norm(time_series-predicted_time_series)/norm(time_series);
    r2_fit_results(rdout_pos(idx),phenc_pos(idx),6)=sum(double(logical(w)));
end

return;