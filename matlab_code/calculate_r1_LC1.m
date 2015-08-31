function [r1_fit_results,tTimeit]=calculate_r1_LC1(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R1 maps using a derivative-free method (the same in NIPALS) on complex-valued data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

tTimeit=0;
startvals=[1.0;1.0;1.0;2.0];
LB=[-5;-5;r1_lb;1];
UB=[5;5;r1_ub;3];

szdata=size(data);
r1_fit_results=zeros([szdata(1:2) 8]);

if isempty(varargin{1})
    temp=ones(szdata(1:2));
    [rdout_pos,phenc_pos]=find(temp);
    clear temp
else
    mask=varargin{1};
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

options=optimset('fminsearch');
options=optimset(options,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off');

for idx=1:length(phenc_pos)
    
    time_series=(squeeze(data(rdout_pos(idx),phenc_pos(idx),:)));
    
    if weight_flag
        w=squeeze(weights(rdout_pos(idx),phenc_pos(idx),:));
    end
    
    if noise_thresh_flag
        w(time_series<noise_thresh)=0;
    end
    
    fit_numvals=sum(double(logical(w)));
    
    if fit_numvals>2
        if idx==1 %benchmarking
            f=@()r1_LC1(time_series,w,ti,startvals,LB,UB,options);
            tTimeit=timeit(f);
        end
        
        [estimates,eflag,fval]=r1_LC1(time_series,w,ti,startvals,LB,UB,options);
        temp=(1-(estimates(4)*exp(-ti(:)*(estimates(4)-1)*estimates(3))));
        
        ILP=estimates(1)+1i*estimates(2);
        predicted=ILP*temp;
        error=norm(time_series-predicted)/norm(time_series);
        
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[estimates(3:4);angle(ILP);abs(ILP);error;eflag;fval;fit_numvals];
    else
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
end

r1_fit_results(:,:,4)=squeeze(r1_fit_results(:,:,4)).*scale_factors;
r1_fit_results=single(r1_fit_results);

%% Function r1_LC1
    function [estimates,eflag,fval]=r1_LC1(data,w,ti,startvals,LB,UB,options)
        ph_weights=(abs(data).^2)./sum(abs(data).^2);
        ph_est=transpose(ph_weights(:))*angle(data(:));
        
        startvals(1)=cos(ph_est);
        startvals(2)=sin(ph_est);
        
        [estimates,eflag,fval]=fminsearchbnd(@(params)obj_r1_LC1(params,ti,data,w),startvals,LB,UB,options);
    end

%% Function obj_r1_LC1
    function obj=obj_r1_LC1(params,time_pts,time_series,w)
        value=(params(1)+1i*params(2))*(1-(params(4)*exp(-time_pts(:)*(params(4)-1)*params(3))));
        obj=norm(w(:).*(value-time_series));
    end
end