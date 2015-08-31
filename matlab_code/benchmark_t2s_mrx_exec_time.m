function [exec_time,t2s_estimate]=benchmark_t2s_mrx_exec_time(data,te,algo)
%benchmark_t2s_mrx_exec_time: Function to benchmark execution times related
%to T2 (T2*) estimation from simulated Multi-rx data.
% Inputs:
% data: A 3d matrix of values, where each row represents a single T2 (T2*)
% data acquired at different TE time points (2nd dimension). The 3rd
% dimension represents various coils.
% te: TE time points (s)
% algo: A number indicating the algorithm to use for execution
% 1:LC 2:VC (see below for details)
% Outputs:
% exec_time: Execution time as calculated by timeit function
% t2_est: Estimated T2 (T2*)

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

% Some initialization for the VARPRO algorithm
global Ind
Ind1=repmat(1:16,[2 1]);
Ind2=repmat(1:2,[1 16]);
Ind=[transpose(Ind1(:));transpose(Ind2(:))];

% Set tolerances for different optimizers
options_lsq=optimset('lsqnonlin');
options_lsq=optimset(options_lsq,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','DerivativeCheck','off','Jacobian','on');

data=squeeze(data);
te=te(:);

% Set lower and upper bounds for the parameters
max_t2star=200e-3;
NLLB=[1/max_t2star;-200];
NLUB=[Inf;200];
NLPstart=[1/max_t2star;0];

ydata=[real(data);imag(data)];
ydata=double(ydata(:));

switch algo
    case 1 % LC
        LPstart=[real(data(1,:));imag(data(1,:))];
        LPstart=LPstart(:);
        LLB=-5*(abs(LPstart));
        LUB=5*(abs(LPstart));
        
        start_vals=[NLPstart(:);LPstart(:)];start_vals=double(start_vals);
        LB=[NLLB(:);LLB(:)];LB=double(LB);
        UB=[NLUB(:);LUB(:)];UB=double(UB);
        
        f=@()lsqnonlin(@(params)t2s_mtx_lsqnonlin(params,te,ydata),start_vals,LB,UB,options_lsq);
        exec_time=timeit(f);
        opt_params=lsqnonlin(@(params)t2s_mtx_lsqnonlin(params,te,ydata),start_vals,LB,UB,options_lsq);
        r2s_estimate=opt_params(1);

    case 2 % VC
        f=@()varpro_basic(ydata, ones(size(ydata)), NLPstart, 16, @(params)t2s_mtx_varpro(params,te), NLLB, NLUB, options_lsq);
        exec_time=timeit(f);
        opt_params=varpro_basic(ydata, ones(size(ydata)), NLPstart, 16, @(params)t2s_mtx_varpro(params,te), NLLB, NLUB, options_lsq);
        r2s_estimate=opt_params(1);
        
    otherwise
        disp('Unknown function')
end

t2s_estimate=1/r2s_estimate;

%% Function t2s_mtx_lsqnonlin
    function [obj,Jmtx]=t2s_mtx_lsqnonlin(params,tes,ydata)
        n=length(tes);
        amps=params(3:end);
        amps=reshape(amps,[2 8]);
        
        vec1=cos(2*pi*params(2)*tes(:));
        vec2=sin(2*pi*params(2)*tes(:));
        
        exp_decay=exp(-tes(:)*params(1));
        
        vec_mtx=[vec1 vec2].*repmat(exp_decay,[1 2]);
        Re_S=vec_mtx*[amps(1,:);-amps(2,:)];
        Im_S=vec_mtx*[amps(2,:);amps(1,:)];
        
        obj=[Re_S;Im_S];
        obj=(obj(:)-ydata);
        
        Ds_by_DR2=[Re_S.*repmat(-tes(:),[1 8]);Im_S.*repmat(-tes(:),[1 8])];
        Ds_by_DR2=Ds_by_DR2(:);
        
        Ds_by_Df=[Im_S.*repmat(-2*pi*tes(:),[1 8]);Re_S.*repmat(2*pi*tes(:),[1 8])];
        Ds_by_Df=Ds_by_Df(:);
        
        Ds_by_DA1=zeros(2*n,1);
        Ds_by_DA1(1:n)=exp_decay.*vec1;
        Ds_by_DA1(n+1:end)=exp_decay.*vec2;
        
        Ds_by_DB1=zeros(2*n,1);
        Ds_by_DB1(1:n)=-exp_decay.*vec2;
        Ds_by_DB1(n+1:end)=exp_decay.*vec1;
        
        temp_vec=[Ds_by_DA1 Ds_by_DB1];
        Jmtx=[Ds_by_DR2 Ds_by_Df blkdiag(temp_vec,temp_vec,temp_vec,temp_vec,temp_vec,temp_vec,temp_vec,temp_vec)];
    end

%% Function t2s_mtx_varpro
    function [Phi,dPhi,Ind]=t2s_mtx_varpro(params,te)
        global Ind
        
        temp_Phi=exp(-te(:)*params(1)).*exp(1i*(2*pi*params(2)*te(:)));
        temp_dPhi=[-temp_Phi(:).*te(:) temp_Phi(:).*(1i*2*pi*te(:))];
        
        temp_Phi=[real(temp_Phi) -imag(temp_Phi);imag(temp_Phi) real(temp_Phi)];
        temp_dPhi=[real(temp_dPhi) -imag(temp_dPhi);imag(temp_dPhi) real(temp_dPhi)];
        
        Phi=(blkdiag(temp_Phi,temp_Phi,temp_Phi,temp_Phi,temp_Phi,temp_Phi,temp_Phi,temp_Phi));
        dPhi=(blkdiag(temp_dPhi,temp_dPhi,temp_dPhi,temp_dPhi,temp_dPhi,temp_dPhi,temp_dPhi,temp_dPhi));
    end
end