function [r2_fit_results,tTimeit,fit_errors]=calculate_r2_LP(data,te,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R2 maps using linear prediction on complex-valued
% data using Steiglitz-McBride iteration
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

delta=mean(diff(te));
szdata=size(data);
r2_fit_results=zeros([szdata(1:2) 8]);
fit_errors=zeros(szdata);
tTimeit=0;

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

%Normalize data for consistent fits
scale_factors=max(abs(data),[],3);
data=data./repmat(scale_factors,[1 1 size(data,3)]);

for idx=1:length(phenc_pos)
    time_series=(squeeze(data(rdout_pos(idx),phenc_pos(idx),:)));
    
    if weight_flag
        w=squeeze(weights(rdout_pos(idx),phenc_pos(idx),:));
    end
    
    if noise_thresh_flag
        w(abs(time_series)<noise_thresh)=0;
    end
    
    fit_numvals=sum(double(logical(w)));
    
    if fit_numvals>3
        if idx==1 %benchmarking
            f=@()stmcb(time_series,0,1);
            tTimeit=timeit(f);
        end
        
        [b,a] = stmcb(time_series,0,1);
        f_est=mod(-imag(log(a(2)))/(2*pi*delta),0.5/delta);
        r2s_est=-(real(log(a(2))))/delta;
        c=b/a(2);
        
        predicted=c*exp(1i*(-2*pi*f_est*te + phase(c))).*exp(-te*r2s_est);
        fit_errors(rdout_pos(idx),phenc_pos(idx),:)=time_series(:)-predicted(:);
        
        error=norm(squeeze(fit_errors(rdout_pos(idx),phenc_pos(idx),:)))/norm(time_series);
        
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[r2s_est;angle(c);f_est;abs(c);error;0;0;fit_numvals];
    else
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
end

r2_fit_results(:,:,4)=squeeze(r2_fit_results(:,:,4)).*scale_factors;
r2_fit_results=single(r2_fit_results);

return;