function [r2_fit_results,tTimeit,fit_errors]=calculate_r2_VM(data,te,t2_upper_threshold,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R2 maps using VARPRO on magnitude data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

NLPstart=10;
INLB=1/t2_upper_threshold;
INUB=Inf;
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
options=optimset(options,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','DerivativeCheck','off','Jacobian','on');

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
        
        if idx==1 % One type of benchmarking
            f=@()varpro_basic(temp_data(:), temp_w(:), NLPstart, 1, @(params)r2_varpro_abs(params,temp_te(:)), INLB, INUB, options);
            tTimeit=timeit(f);
        end
        
        [INLP, ILP, ~, wresid_norm, eflag] = varpro_basic(temp_data(:), temp_w(:), NLPstart, 1, @(params)r2_varpro_abs(params,temp_te(:)), INLB, INUB, options);
        
        predicted_time_series=ILP*exp(-te(:)*INLP);
        
        fit_errors(rdout_pos(idx),phenc_pos(idx),:)=time_series(:)-predicted_time_series(:);
        error=norm(squeeze(fit_errors(rdout_pos(idx),phenc_pos(idx),:)))/norm(time_series);
        
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[INLP;ILP;error;eflag;wresid_norm;fit_numvals];
    else
        tTimeit=0;
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
end
r2_fit_results(:,:,2)=squeeze(r2_fit_results(:,:,2)).*scale_factors;
r2_fit_results=single(r2_fit_results);


%% Function r2_varpro_abs
    function [Phi,dPhi,Ind]=r2_varpro_abs(params,te)
        % See varpro_basic for structure of this function
        Phi=exp(-te(:)*params(1));
        Ind=[1;1];
        dPhi=-Phi(:).*te(:);
    end
end