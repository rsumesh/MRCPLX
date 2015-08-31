function [exec_time,t1_est]=benchmark_mrx_t1_exec_time(data,ti,algo)
%benchmark_mrx_t1_exec_time: Function to benchmark execution times related
%to T1 estimation from simulated Multi-rx data.
% Inputs:
% data: A 3d matrix of values, where each row represents a single T1
% data acquired at different TI time points (2nd dimension). The 3rd
% dimension represents various coils.
% ti: TI (inversion) time points (s)
% algo: A number indicating the algorithm to use for execution
% 1:LC 2:VC (see below for details)
% Outputs:
% exec_time: Execution time as calculated by timeit function
% t1_est: Estimated T1

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY


% Some initialization for the VARPRO algorithm
global Jmtx Phi dPhi Ind

n=length(ti);
ncoils=size(squeeze(data),2);

Jmtx=zeros(2*numel(data),2+2*ncoils);

sub1_Jmtx_A=1:2*n*ncoils;
sub2_Jmtx_A=repmat(3:2:(2+2*ncoils),[2*n 1]);
sub2_Jmtx_B=repmat(4:2:(2+2*ncoils),[2*n 1]);

sub1_Jmtx_A=reshape(sub1_Jmtx_A(:),[2*n ncoils]);
sub2_Jmtx_A=reshape(sub2_Jmtx_A(:),[2*n ncoils]);
sub2_Jmtx_B=reshape(sub2_Jmtx_B(:),[2*n ncoils]);

ind_Jmtx_A=sub2ind([size(Jmtx,1) size(Jmtx,2)],sub1_Jmtx_A,sub2_Jmtx_A);
ind_Jmtx_B=sub2ind([size(Jmtx,1) size(Jmtx,2)],sub1_Jmtx_A,sub2_Jmtx_B);

Ind1=repmat(1:2*ncoils,[2 1]);
Ind2=repmat(1:2,[1 2*ncoils]);
Ind=[transpose(Ind1(:));transpose(Ind2(:))];

Phi=zeros(2*numel(data),2*ncoils);
dPhi=zeros(2*numel(data),4*ncoils);

sub1_Phi=1:n*2*ncoils;
sub2_Phi=repmat(1:2*ncoils,[n 1]);

sub2_dPhi_D=repmat([1:2:4*ncoils],[n 1]);
sub2_dPhi_C=repmat([2:2:4*ncoils],[n 1]);

sub1_Phi=reshape(sub1_Phi,[n 2*ncoils]);
sub2_Phi=reshape(sub2_Phi,[n 2*ncoils]);
sub2_dPhi_D=reshape(sub2_dPhi_D,[n 2*ncoils]);
sub2_dPhi_C=reshape(sub2_dPhi_C,[n 2*ncoils]);

ind_Phi=sub2ind([size(Phi,1) size(Phi,2)],sub1_Phi,sub2_Phi);
ind_dPhi_D=sub2ind([size(dPhi,1) size(dPhi,2)],sub1_Phi,sub2_dPhi_D);
ind_dPhi_C=sub2ind([size(dPhi,1) size(dPhi,2)],sub1_Phi,sub2_dPhi_C);

% Set tolerances for different optimizers
options_lsq=optimset('lsqnonlin');
options_lsq=optimset(options_lsq,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e4,'Display','off','DerivativeCheck','off','Jacobian','on');

% Set lower and upper bounds for the parameters
r1_lb=0.1;
r1_ub=10;
NLPstart=[1.0;2.0];
NLLB=[r1_lb;1.0];
NLUB=[r1_ub;3.0];

ti=ti(:);

% Scale the data
data=squeeze(data)/max(abs(data(:)));

switch algo
    case 1 % LC-CPLX MRX
        f=@()r1_LC_mrx(data,ti,NLPstart,NLLB,NLUB,options_lsq);
        exec_time=timeit(f);
        r1_est=r1_LC_mrx(data,ti,NLPstart,NLLB,NLUB,options_lsq);
        
    case 2 % VC-CPLX MRX
        ydata=[real(data(:));imag(data(:))];
        f=@()varpro_basic(ydata, ones(size(ydata)), NLPstart, 2*ncoils, @(params)r1_VC_mrx(params,ti), NLLB, NLUB, options_lsq);
        exec_time=timeit(f);
        INLP=varpro_basic(ydata, ones(size(ydata)), NLPstart, 2*ncoils, @(params)r1_VC_mrx(params,ti), NLLB, NLUB, options_lsq);
        temp_C=INLP(2);
        temp_D=INLP(1);
        r1_est=temp_D/(temp_C -1);
        
    otherwise
        disp('Unknown function')
end

t1_est=1/r1_est;

%% Function r1_LC
    function [r1_est,estimates,fval,eflag]= r1_LC_mrx(data,ti,startvals,LB,UB,options)
        ph_weights=(abs(data).^2)./repmat(sum(abs(data).^2),[size(data,1) 1]);
        ph_est=diag(transpose(ph_weights)*angle(data));
        
        startvals_real=cos(ph_est);
        startvals_imag=sin(ph_est);
        temp_startvals=[startvals_real;startvals_imag];
        temp_startvals=[startvals;temp_startvals(:)];
        
        LB1=[LB;-5*ones(2*ncoils,1)];
        UB1=[UB;5*ones(2*ncoils,1)];
        
        ydata1=[real(data);imag(data)];
        ydata1=ydata1(:);
        
        [estimates,fval,~,eflag]=lsqnonlin(@(params)obj_r1_LC_mrx(params,ti,ydata1),temp_startvals,LB1,UB1,options);
        r1_est=estimates(1)/(estimates(2)-1);
    end

%% Function obj_r1_LC
    function [obj,Jmtx]=obj_r1_LC_mrx(params,ti,ydata)
        
        global Jmtx

        D=params(1);
        C=params(2);
        
        amps=params(3:end);
        amps=reshape(amps,[2 ncoils]);
        amps=amps(1,:)+1i*amps(2,:);
        
        exp_part=exp(-ti(:)*D);
        temp_Phi=(1 - C*exp_part);
        
        Jvec_D=(C*exp_part.*ti(:))*amps;
        Jvec_D=[real(Jvec_D);imag(Jvec_D)];
        Jvec_D=Jvec_D(:);
        
        Jvec_C=-exp_part*amps;
        Jvec_C=[real(Jvec_C);imag(Jvec_C)];
        Jvec_C=Jvec_C(:);
        
        dS_by_DA1=[temp_Phi;zeros(length(ti),1)];
        dS_by_DB1=[zeros(length(ti),1);temp_Phi];

        Jmtx(:,1)=Jvec_D;
        Jmtx(:,2)=Jvec_C;
        
        Jmtx(ind_Jmtx_A(:))=repmat(dS_by_DA1,[size(ind_Jmtx_A,2) 1]);
        Jmtx(ind_Jmtx_B(:))=repmat(dS_by_DB1,[size(ind_Jmtx_A,2) 1]);
                
        temp_mtx=temp_Phi*amps;
        temp_mtx=[real(temp_mtx);imag(temp_mtx)];
        obj=(temp_mtx(:)-ydata);
    end

%% Function r1_VC
    function [Phi,dPhi,Ind]=r1_VC_mrx(params,ti)
        global Phi dPhi Ind
        
        D=params(1);
        C=params(2);
        
        exp_part=exp(-ti(:)*D);
        temp_Phi=(1 - C*exp_part);
       
        temp_dPhi_D=(C*exp_part.*ti(:));
        temp_dPhi_C=-exp_part;
        
        Phi(ind_Phi(:))=repmat(temp_Phi,[size(ind_Phi,2) 1]);
        dPhi(ind_dPhi_D(:))=repmat(temp_dPhi_D,[size(ind_Phi,2) 1]);
        dPhi(ind_dPhi_C(:))=repmat(temp_dPhi_C,[size(ind_Phi,2) 1]);
    end
end