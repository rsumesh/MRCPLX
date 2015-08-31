function [r2_fit_results,tTimeit,fit_errors]=calculate_r2_LQ(data,te,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R2 maps using Levenberg-Marquardt method on magnitude data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

start_vals=[1;10];
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

%Normalize data for consistent fits
scale_factors=max(abs(data),[],3);
data=data./repmat(scale_factors,[1 1 size(data,3)]);

options=optimset('lsqnonlin');
options=optimset(options,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','algorithm','levenberg-marquardt','DerivativeCheck','off','Jacobian','on');

for idx=1:length(phenc_pos)
    time_series=abs(squeeze(data(rdout_pos(idx),phenc_pos(idx),:)));
    
    if weight_flag
        w=squeeze(weights(rdout_pos(idx),phenc_pos(idx),:));
    end
    if noise_thresh_flag
        w(abs(time_series)<noise_thresh)=0;
    end
    
    fit_numvals=sum(double(logical(w)));
    
    if fit_numvals>3
        
        
        temp=find(w==0);
        if ~isempty(temp) && temp(1)>2
            temp_data=time_series(1:(temp(1)-1));
            temp_w=w(1:(temp(1)-1));
            temp_te=te(1:(temp(1)-1));
        else
            temp_data=time_series;
            temp_w=w;
            temp_te=te;
        end
        
        if idx==1 %benchmarking
            f=@()lsqnonlin(@(params)r2_LMQ(params,temp_te(:),temp_data(:),temp_w),start_vals,[],[],options);
            tTimeit=timeit(f);
        end
        
        [opt_params,fval,~,eflag]=lsqnonlin(@(params)r2_LMQ(params,temp_te(:),temp_data(:),temp_w),start_vals,[],[],options);
        
        predicted=opt_params(1)*exp(-te(:)*opt_params(2));
        
        fit_errors(rdout_pos(idx),phenc_pos(idx),:)=time_series(:)-predicted(:);
        
        error=norm(squeeze(fit_errors(rdout_pos(idx),phenc_pos(idx),:)))/norm(time_series);
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[opt_params(2);opt_params(1);error;eflag;abs(fval);fit_numvals];
    else
        tTimeit=0;
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
    
end
r2_fit_results(:,:,4)=squeeze(r2_fit_results(:,:,4)).*scale_factors;
r2_fit_results=single(r2_fit_results);

%% Function r2_LMQ
    function [obj,J]=r2_LMQ(params,te,time_series,w)
        
        Phi=exp(-te(:)*params(2));
        obj=(w(:).*(params(1)*Phi-time_series));% Objective
        J=[Phi -params(1)*Phi.*te(:)];% Jacobian
    end
end