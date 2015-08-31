function [r1_fit_results,tTimeit]=calculate_r1_LQ(data,ti,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R1 maps using the Levenberg-Marquardt method
% Since derivatives do not exist for magnitude T1-data, finite difference
% approximations have been used.
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

tTimeit=0;
startvals=[1.0;1.0;2.0];

szdata=size(data);
r1_fit_results=zeros([szdata(1:2) 7]);

if isempty(varargin{1})
    temp=ones(szdata(1:2));
    [rdout_pos,phenc_pos]=find(temp);
    clear temp
else
    mask=varargin{1,1};
    [rdout_pos,phenc_pos]=find(mask);
end

if weight_flag
    weights=abs(data);
else
    w=ones(size(data,3),1);
end

%Normalize data for consistent fits
scale_factors=max(abs(data),[],3);
data=data./repmat(scale_factors,[1 1 size(data,3)]);

options=optimset('lsqnonlin');
options=optimset(options,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','algorithm','levenberg-marquardt');

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
        
        if idx==1 %benchmarking
            f=@()lsqnonlin(@(params)r1_LMQ(params,ti,time_series,w),startvals,[],[],options);
            tTimeit=timeit(f);
        end
        
        [estimates,fval,~,eflag]=lsqnonlin(@(params)r1_LMQ(params,ti,time_series,w),startvals,[],[],options);
        
        predicted=abs(estimates(1)*(1-(estimates(3)*exp(-ti(:)*(estimates(3)-1)*estimates(2)))));
        error=norm(time_series-predicted)/norm(time_series);
        
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[estimates(2);estimates(3);estimates(1);error;eflag;fval;fit_numvals];
    else
        tTimeit=Inf;
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
end

r1_fit_results(:,:,3)=squeeze(r1_fit_results(:,:,3)).*scale_factors;
r1_fit_results=single(r1_fit_results);

%% Function r1_LMQ
    function obj=r1_LMQ(params,ti,ydata,w)
        R1=params(2);
        B=params(3);
        exp_part=exp(-ti(:)*(B-1)*R1);
        
        obj=w(:).*(ydata(:) - abs(params(1)*(1- B*exp_part)));
    end
end