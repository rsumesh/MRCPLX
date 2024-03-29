function [r1_fit_results,tTimeit]=calculate_r1_VC(data,ti,r1_lb,r1_ub,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R2 maps using VARPRO on complex-valued data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

tTimeit=0;
NLPstart=1.0;
INLB=r1_lb;
INUB=r1_ub;

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
    weights=cat(3,abs(real(data)),abs(imag(data)));
else
    w=ones(2*size(data,3),1);
end

%Normalize data for consistent fits
scale_factors=max(abs(data),[],3);
data=data./repmat(scale_factors,[1 1 size(data,3)]);

options=optimset('lsqnonlin');
options = optimset(options,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','DerivativeCheck','off','Jacobian','on');

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
        ydata=[real(time_series(:));imag(time_series(:))];
        
        if idx==1 %benchmarking
            f=@()varpro_basic(ydata, w, NLPstart, 4, @(params)r1_VC_type3(params,ti), INLB, INUB, options);
            tTimeit=timeit(f);
        end
        
        [INLP, ILP, ~, wresid_norm, exit_flag] = varpro_basic(ydata, w, NLPstart, 4, @(params)r1_VC_type3(params,ti), INLB, INUB, options);
        
        c1=ILP(1)+1i*ILP(2);
        c2=ILP(3)+1i*ILP(4);
        
        temp_r1=INLP/(abs(c2/c1)-1);
        
        predicted=c1 + c2*exp(-ti(:)*INLP);
        
        error=norm(time_series(:)-predicted(:))/norm(time_series);
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[temp_r1;abs(c2/c1);angle((ILP(1)+(1i*ILP(2))));abs((ILP(1)+(1i*ILP(2))));error;exit_flag;wresid_norm;fit_numvals];
    else
        r1_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
end

r1_fit_results(:,:,4)=squeeze(r1_fit_results(:,:,4)).*scale_factors;
r1_fit_results=single(r1_fit_results);

%% Function r1_VC_type3
    function [Phi,dPhi,Ind]=r1_VC_type3(params,ti)
        exp_part=exp(-ti(:)*params(1));
        
        zz=zeros(size(exp_part,1),1);
        zz1=zz+1;
        
        Phi=[zz1 zz exp_part(:) zz;zz zz1 zz exp_part(:)];
        
        Ind=[3 4;1 1];
        temp=-ti(:).*exp_part;
        dPhi=[temp zz;zz temp];
    end

end