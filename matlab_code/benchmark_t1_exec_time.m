function [exec_time,t1_est]=benchmark_t1_exec_time(data,ti,algo)
%benchmark_t1_exec_time: Function to benchmark processing times for various
%T1 estimation algorithms
% Inputs:
% data: A matrix of values, where each row represents a single T1
% data acquired at different TI time points.
% ti: TI time points (s)
% algo: A number indicating the algorithm to use for execution
% 1: LQ, 2: NM, 3: LC, 4: NC, 5:VC (see below for details)
% Outputs:
% exec_time: Execution time as calculated by timeit function
% t1_est: Estimated T1

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

% Set tolerances for different optimizers
options_lsq=optimset('lsqnonlin');
options_lsq=optimset(options_lsq,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off');
options_lsq_lmq=optimset(options_lsq,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','algorithm','levenberg-marquardt');
options_lsq=optimset(options_lsq,'DerivativeCheck','off','Jacobian','on');
options_fmins=optimset('fminsearch');
options_fmins=optimset(options_fmins,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','Algorithm','interior-point');

% Set lower and upper bounds for the parameters
r1_lb=0.1;
r1_ub=10;
w=ones(size(data));
startvals=[1.0;1.0;2.0];
NLPstart=[1.0;2.0];
NLLB=[r1_lb;1.0];
NLUB=[r1_ub;3.0];

cplx_startvals=[1.0;1.0;1.0;2.0];
cplx_LB=[-5;-5;r1_lb;1];
cplx_UB=[5;5;r1_ub;3];

data=data(:);
ti=ti(:);


data=data/max(abs(data));%Scale data

switch algo
    case 1 % LQ-MAG Levenberg-Marquardt fit
        f=@()lsqnonlin(@(params)r1_LQ(params,ti,abs(data),w),startvals,[],[],options_lsq_lmq);
        exec_time=timeit(f);
        estimates=lsqnonlin(@(params)r1_LQ(params,ti,abs(data),w),startvals,[],[],options_lsq_lmq);
        r1_est=estimates(2);
        
    case 2 % NM-MAG NIPALS Gaussian Assumption
        f=@()fminspleas({@(params,ti)abs(1-(params(2)*exp(-ti(:)*(params(2)-1)*params(1))))},NLPstart,ti(:),abs(data),NLLB,NLUB,w,options_fmins);
        exec_time=timeit(f);
        INLP=fminspleas({@(params,ti)abs(1-(params(2)*exp(-ti(:)*(params(2)-1)*params(1))))},NLPstart,ti(:),abs(data),NLLB,NLUB,w,options_fmins);
        r1_est=INLP(1);
        
    case 3 % LC-CPLX % Derivative-based Least-squares
        temp_w=ones(2*length(w),1);
        f=@()r1_LC(data,temp_w,ti,cplx_startvals,cplx_LB,cplx_UB,options_lsq);
        exec_time=timeit(f);
        r1_est=r1_LC(data,temp_w,ti,cplx_startvals,cplx_LB,cplx_UB,options_lsq);
        
    case 4 % NC-CPLX % NIPALS on complex-valued data
        f=@()fminspleas({@(params,ti)(1-(params(2)*exp(-ti(:)*(params(2)-1)*params(1))))},NLPstart,ti,data,NLLB,NLUB,w,options_fmins);
        exec_time=timeit(f);
        INLP=fminspleas({@(params,ti)(1-(params(2)*exp(-ti(:)*(params(2)-1)*params(1))))},NLPstart,ti,data,NLLB,NLUB,w,options_fmins);
        r1_est=INLP(1);
        
    case 5 % VC-CPLX % VARPRO on complex-valued data
        temp_w=ones(2*length(w),1);
        ydata=[real(data(:));imag(data(:))];
        f=@()varpro_basic(ydata, temp_w, NLPstart(1), 4, @(params)r1_VC(params,ti), r1_lb, r1_ub, options_lsq);
        exec_time=timeit(f);
        [INLP, ILP]=varpro_basic(ydata, temp_w, NLPstart(1), 4, @(params)r1_VC(params,ti), r1_lb, r1_ub, options_lsq);
        
        c1=ILP(1)+1i*ILP(2);
        c2=ILP(3)+1i*ILP(4);
        
        r1_est=INLP/(abs(c2/c1)-1);
        
    otherwise
        disp('Unknown function')
end

t1_est=1/r1_est;

%% Function r1_LQ
    function obj=r1_LQ(params,ti,ydata,w)
        R1=params(2);
        B=params(3);
        exp_part=exp(-ti(:)*(B-1)*R1);
        
        obj=w(:).*(ydata(:) - abs(params(1)*(1- B*exp_part)));
    end

%% Function r1_LC
    function [r1_est,estimates,fval,eflag]= r1_LC(data,w,ti,startvals,LB,UB,options)
        ph_weights=(abs(data).^2)./sum(abs(data).^2);
        ph_est=transpose(ph_weights(:))*angle(data(:));
        
        startvals(1)=cos(ph_est);
        startvals(2)=sin(ph_est);
        
        ydata1=[real(data(:));imag(data(:))];
        
        [estimates,fval,~,eflag]=lsqnonlin(@(params)obj_r1_LC(params,ti,ydata1,w),startvals,LB,UB,options);
        r1_est=estimates(3)/(estimates(4)-1);
    end

%% Function obj_r1_LC
    function [obj,Jmtx]=obj_r1_LC(params,ti,ydata,w)
        D=params(3);
        C=params(4);
        
        amp=params(1) + 1i*params(2);
        exp_part=exp(-ti(:)*D);
        Phi=(1 - C*exp_part);
        
        Jmtx=zeros(size(ti,1),4);
        Jmtx(:,1)=Phi;
        Jmtx(:,2)=1i*Phi;
        Jmtx(:,3)=amp*(Phi -1).*(-ti(:));
        Jmtx(:,4)=-amp*exp_part;
        
        Jmtx=[real(Jmtx);imag(Jmtx)];
        obj=w(:).*([params(1)*Phi;params(2)*Phi]-ydata);
    end

%% Function r1_VC
    function [Phi,dPhi,Ind]=r1_VC(params,ti)
        exp_part=exp(-ti(:)*params(1));
        
        zz=zeros(size(exp_part,1),1);
        zz1=zz+1;
        
        Phi=[zz1 zz exp_part(:) zz;zz zz1 zz exp_part(:)];
        
        Ind=[3 4;1 1];
        temp=-ti(:).*exp_part;
        dPhi=[temp zz;zz temp];
    end
end