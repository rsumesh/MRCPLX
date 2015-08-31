function [r2_fit_results,tTimeit,fit_errors]=calculate_r2_LC(data,te,t2_upper_threshold,max_f_deviation,weight_flag,noise_thresh_flag,noise_thresh,varargin)

% Function to calculate R2 maps using a derivative-based regular
% least-squares method on complex-valued data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

start_vals=[1;0;10;0];
LB=[-5;-5;1/t2_upper_threshold;-max_f_deviation];
UB=[5;5;Inf;max_f_deviation];
tTimeit=0;

szdata=size(data);
r2_fit_results=zeros([szdata(1:2) 8]);
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

H=[ones(length(te),1) te(:)];
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
     
        ph_weights=(abs(time_series).^2)./sum(abs(time_series).^2);
        ph_wmtx=diag(ph_weights);
        ph_params=(ph_wmtx*H)\(ph_wmtx*(unwrap(angle(time_series))));
        
        % Check if start_vals phase is outside -pi to pi, in that case,
        % bring it back to the right range!!
        if ph_params(1) < -pi
            ph_params(1)=2*pi+ph_params(1);
        end
        if ph_params(1) > pi
            ph_params(1)=2*pi-ph_params(1);
        end
        
        start_vals(1)=cos(ph_params(1));
        start_vals(2)=sin(ph_params(1));
        start_vals(4)=ph_params(2)/(2*pi);
        
        if idx==1 %benchmarking
            f=@()lsqnonlin(@(params)r2_LC2(params,te,[real(time_series(:));imag(time_series(:))],[w(:);w(:)]),start_vals,LB,UB,options);
            tTimeit=timeit(f);
        end
        
        [opt_params,fval,~,eflag]=lsqnonlin(@(params)r2_LC2(params,te,[real(time_series(:));imag(time_series(:))],[w(:);w(:)]),start_vals,LB,UB,options);
        
        ILP=opt_params(1)+1i*opt_params(2);
        predicted=ILP*exp(-te(:)*opt_params(3)).*exp(1i*2*pi*opt_params(4)*te(:));
        
        fit_errors(rdout_pos(idx),phenc_pos(idx),:)=time_series(:)-predicted(:);
        
        error=norm(squeeze(fit_errors(rdout_pos(idx),phenc_pos(idx),:)))/norm(time_series);
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[opt_params(3);angle(ILP);opt_params(4);abs(ILP);error;eflag;abs(fval);fit_numvals];
    else
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
end

r2_fit_results(:,:,4)=squeeze(r2_fit_results(:,:,4)).*scale_factors;
r2_fit_results=single(r2_fit_results);

%% Function r2_LC2
    function [obj,J]=r2_LC2(params,te,time_series,w)
        
        Phi=exp(-te(:)*params(3)).*exp(1i*(2*pi*params(4)*te(:)));
        amp=params(1)+1i*params(2);
        
        J=[Phi 1i*Phi amp*Phi.*-te(:) amp*Phi.*(1i*2*pi*te(:))];
        J=[real(J);imag(J)];% Jacobian
        
        magnetization=amp*Phi;
        magnetization=[real(magnetization);imag(magnetization)];
        
        obj=w.*(magnetization-time_series);%Objective
    end
end