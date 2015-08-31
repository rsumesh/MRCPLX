function [exec_time,t2_est]=benchmark_t2_exec_time(data,te,algo)
%benchmark_t2_exec_time: Function to benchmark processing times for various
%T2 (T2*) estimation algorithms
% Inputs:
% data: A matrix of values, where each row represents a single T2 (T2*)
% data acquired at different TE time points.
% te: TE time points (s)
% algo: A number indicating the algorithm to use for execution
% 1:LG, 2:LQ, 3:AR, 4:NM, 5:LC, 6:NC, 7:VC, 8:LP (see below for details)
% Outputs:
% exec_time: Execution time as calculated by timeit function
% t2_est: Estimated T2 (T2*)

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

% Set lower and upper bounds for the parameters
r2_upper_threshold=200;
max_f_deviation=0.5/max(diff(te(:)));

start_vals=[1;10];
lb=[0;0];
ub=[5;r2_upper_threshold];

start_vals_cplx=[1;0;10;0];
lb_cplx=[-5;-5;0;-max_f_deviation];
ub_cplx=[5;5;r2_upper_threshold;max_f_deviation];

data=data(:);
te=te(:);
w=ones(size(data));
N=length(data);

% Set tolerances for different optimizers
options_lsq=optimset('lsqnonlin');
options_lsq_lmq=optimset(options_lsq,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','algorithm','levenberg-marquardt','DerivativeCheck','off','Jacobian','on');
options_lsq=optimset(options_lsq,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','DerivativeCheck','off','Jacobian','on');

options_fmins=optimset('fminsearch');
options_fmins=optimset(options_fmins,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','Algorithm','interior-point');

% Scale the data
data=data/max(abs(data));

switch algo
    case 1 % LG-MAG Linear fit using log operation on magnitude data
        H=[ones(N,1) -te];
        W_mtx=diag(w.^0.5);
        pinvH=pinv(W_mtx*H);
        
        f=@()r2_LG(abs(data),pinvH,W_mtx);
        exec_time=timeit(f);
        
        estimates=r2_LG(abs(data),pinvH,W_mtx);
        r2_est=estimates(2);
        
    case 2 % LQ-MAG Levenberg-Marquardt fit to magnitude data
        f=@()lsqnonlin(@(params)r2_LQ(params,te,abs(data),w),start_vals,[],[],options_lsq_lmq);
        exec_time=timeit(f);
        estimates=lsqnonlin(@(params)r2_LQ(params,te,abs(data),w),start_vals,[],[],options_lsq_lmq);
        r2_est=estimates(2);
        
    case 3 % AR-MAG ARLO fit to magnitude data
        delta=mean(diff(te));
        alpha=(delta/3);
        M=zeros(3,length(te)-2);
        for idx=1:(length(te)-2)
            M(:,idx)=idx:(idx+2);
        end
        f=@()r2_AR(abs(data),M,alpha,te);
        exec_time=timeit(f);
        r2_est=r2_AR(abs(data),M,alpha,te);
        
    case 4 % NM-MAG NIPALS Gaussian Assumption
        f=@()fminspleas({@(param,xdata)exp(-xdata*param)},start_vals(2),te,abs(data),lb(2),ub(2),w,options_fmins);
        exec_time=timeit(f);
        r2_est=fminspleas({@(param,xdata)exp(-xdata*param)},start_vals(2),te,abs(data),lb(2),ub(2),w,options_fmins);
        
        
    case 5 % LC-CPLX Derivative-based Least squares
        H=[ones(length(te),1) te];
        f=@()r2_LC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_lsq);
        exec_time=timeit(f);
        estimates=r2_LC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_lsq);
        r2_est=estimates(3);
        
    case 6 % NC-CPLX NIPALS fit to complex-valued data
        H=[ones(length(te),1) te];
        f=@()r2_NC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_fmins);
        exec_time=timeit(f);
        estimates=r2_NC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_fmins);
        r2_est=estimates(1);
        
    case 7 % VC-CPLX VARPRO fit to complex-valued data
        H=[ones(length(te),1) te];
        f=@()r2_VC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_lsq);
        exec_time=timeit(f);
        estimates=r2_VC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_lsq);
        r2_est=estimates(1);
        
    case 8 % LP-CPLX Linear prediction using Steiglitz-McBride iteration
        f=@()stmcb(data,0,1);
        exec_time=timeit(f);
        [b,a] = stmcb(data,0,1);
        delta=mean(diff(te));
        r2_est=-(real(log(a(2))))/delta;
        
    otherwise
        disp('Unknown Algorithm!');
end

t2_est=1/r2_est;

%% Function LG
    function estimates=r2_LG(data,pinvH,W_mtx)
        log_time_series=W_mtx*log(data);
        estimates=pinvH*log_time_series;
    end

%% Function r2_LQ
    function [obj,J]=r2_LQ(params,te,data,w)
        Phi=exp(-te*params(2));
        obj=(w.*(params(1)*Phi-data));
        J=[Phi -params(1)*Phi.*te];
    end

%% Function r2_AR
    function [r2s_est,A,nonlin_part]=r2_AR(data,M,alpha,te)
        temp1=data(M);
        
        S=alpha*(temp1(1,:)+4*temp1(2,:)+temp1(3,:));
        d=temp1(1,:)-temp1(3,:);
        
        p=sum(S.^2);
        q=sum(S.*d);
        r2s_est=(alpha*p + q)/(p + alpha*q);
        nonlin_part=exp(-te*r2s_est);
        A=nonlin_part\data;
    end


%% Function r2_LC
    function estimates=r2_LC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_lsq)
        ph_weights=(abs(data).^2)./sum(abs(data).^2);
        ph_wmtx=diag(ph_weights);
        ph_params=(ph_wmtx*H)\(ph_wmtx*(unwrap(angle(data))));
        
        % Check if start_vals phase is outside -pi to pi, in that case,
        % bring it back to the right range
        if ph_params(1) < -pi
            ph_params(1)=2*pi+ph_params(1);
        end
        if ph_params(1) > pi
            ph_params(1)=2*pi-ph_params(1);
        end
        
        start_vals_cplx(1)=cos(ph_params(1));
        start_vals_cplx(2)=sin(ph_params(1));
        start_vals_cplx(4)=ph_params(2)/(2*pi);
        
        estimates=lsqnonlin(@(params)obj_r2_LC(params,te,[real(data);imag(data)],[w;w]),start_vals_cplx,lb_cplx,ub_cplx,options_lsq);
    end

%% Function obj_r2_LC
    function [obj,J]=obj_r2_LC(params,te,data,w)
        
        Phi=exp(-te*params(3)).*exp(1i*(2*pi*params(4)*te));
        amp=params(1)+1i*params(2);
        
        J=[Phi 1i*Phi amp*Phi.*-te amp*Phi.*(1i*2*pi*te)];
        J=[real(J);imag(J)];
        
        magnetization=amp*Phi;
        magnetization=[real(magnetization);imag(magnetization)];
        
        obj=w.*(magnetization-data);
    end

%% Function r2_NC
    function estimates=r2_NC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_fmins)
        ph_weights=(abs(data).^2)./sum(abs(data).^2);
        ph_wmtx=diag(ph_weights);
        ph_params=(ph_wmtx*H)\(ph_wmtx*(unwrap(angle(data))));
        
        % Check if start phase is outside -pi to pi, in that case,
        % bring it back to the right range!!
        if ph_params(1) < -pi
            ph_params(1)=2*pi+ph_params(1);
        end
        if ph_params(1) > pi
            ph_params(1)=2*pi-ph_params(1);
        end
        
        start_vals_cplx=start_vals_cplx(3:4);
        lb_cplx=lb_cplx(3:4);
        ub_cplx=ub_cplx(3:4);
        
        start_vals_cplx(2)=ph_params(2)/(2*pi);
        estimates=fminspleas({@(params,te)exp(-te*params(1)).*exp(1i*2*pi*params(2)*te)},start_vals_cplx,te,data,lb_cplx,ub_cplx,w,options_fmins);
    end

%% Function r2_VC
    function estimates=r2_VC(data,te,H,start_vals_cplx,lb_cplx,ub_cplx,w,options_lsq)
        ph_weights=(abs(data).^2)./sum(abs(data).^2);
        ph_wmtx=diag(ph_weights);
        ph_params=(ph_wmtx*H)\(ph_wmtx*(unwrap(angle(data))));
        
        % Check if start phase is outside -pi to pi, in that case,
        % bring it back to the right range!!
        if ph_params(1) < -pi
            ph_params(1)=2*pi+ph_params(1);
        end
        if ph_params(1) > pi
            ph_params(1)=2*pi-ph_params(1);
        end
        
        start_vals_cplx=start_vals_cplx(3:4);
        lb_cplx=lb_cplx(3:4);
        ub_cplx=ub_cplx(3:4);
        
        start_vals_cplx(2)=ph_params(2)/(2*pi);
        estimates = varpro_basic([real(data);imag(data)], [w;w], start_vals_cplx, 2, @(params)obj_r2_VC(params,te), lb_cplx, ub_cplx, options_lsq);
    end

%% Function obj_r2_VC
    function [Phi,dPhi,Ind]=obj_r2_VC(params,te)
        
        Phi=exp(-te(:)*params(1)).*exp(1i*(2*pi*params(2)*te(:)));
        
        Ind=[1 1 2 2;1 2 1 2];
        
        dPhi=[-Phi(:).*te(:) Phi(:).*(1i*2*pi*te(:))];
        
        Phi=[real(Phi) -imag(Phi);imag(Phi) real(Phi)];
        dPhi=[real(dPhi) -imag(dPhi);imag(dPhi) real(dPhi)];
    end
end