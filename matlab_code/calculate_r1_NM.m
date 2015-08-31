function [r1_fit_results,tTimeit]=calculate_r1_NM(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R1 maps using NIPALS on magnitude data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

tTimeit=0;
NLPstart=[1.0;2.0];
INLB=[r1_lb;1.0];
INUB=[r1_ub;3.0];

szdata=size(data);
r1_fit_results=zeros([szdata(1:2) 7]);

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
            f=@()fminspleas({@(params,ti)abs(1-(params(2)*exp(-ti(:)*(params(2)-1)*params(1))))},NLPstart,ti(:),time_series,INLB,INUB,w,options);
            tTimeit=timeit(f);
        end
        
        [INLP,ILP,eflag,fval]=fminspleas({@(params,ti)abs(1-(params(2)*exp(-ti(:)*(params(2)-1)*params(1))))},NLPstart,ti(:),time_series,INLB,INUB,w,options);
        predicted_time_series=ILP*abs(1-(INLP(2)*exp(-ti(:)*(INLP(2)-1)*INLP(1))));
        
        error=norm(time_series-predicted_time_series)/norm(time_series);
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[INLP(:);ILP;error;eflag;fval;fit_numvals];
    else
        tTimeit=0;
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
end
r1_fit_results(:,:,3)=squeeze(r1_fit_results(:,:,3)).*scale_factors;
r1_fit_results=single(r1_fit_results);

return;