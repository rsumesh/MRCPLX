function [r2_fit_results,tTimeit,fit_errors]=calculate_r2_AR(data,te,weight_flag,noise_thresh_flag,noise_thresh,varargin)
% Function to estimate R2 using ARLO method on magnitude data
% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY

delta=mean(diff(te));
alpha=(delta/3);
szdata=size(data);
r2_fit_results=zeros([szdata(1:2) 6]);
fit_errors=zeros(szdata);
tTimeit=0;

M=zeros(3,length(te)-2);
for idx=1:(length(te)-2)
    M(:,idx)=idx:(idx+2);
end

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

for idx=1:length(phenc_pos)
    time_series=abs(squeeze(data(rdout_pos(idx),phenc_pos(idx),:)));
    
    if weight_flag
        w=squeeze(weights(rdout_pos(idx),phenc_pos(idx),:));
    end
    if noise_thresh_flag
        w(abs(time_series)<noise_thresh)=0;
    end
    
    fit_numvals=sum(double(logical(w)));
    
    if fit_numvals>3
        temp=find(w==0);
        if ~isempty(temp) && temp(1)>2
            temp_data=time_series(1:(temp(1)-1));
            temp_M=M(:,1:(length(temp_data)-2));
            temp_te=te(1:length(temp_data));
        else
            temp_data=time_series;
            temp_M=M;
            temp_te=te;
        end
        
        if idx==1 %benchmarking
            f=@()r2_ARL(temp_data,temp_M,alpha,temp_te);
            tTimeit=timeit(f);
        end
        
        [r2s_est,A,nonlin_part]=r2_ARL(temp_data,temp_M,alpha,temp_te);
        predicted=A*nonlin_part;
        
        fit_errors(rdout_pos(idx),phenc_pos(idx),1:length(temp_data))=temp_data(:)-predicted(:);
        error=norm(squeeze(fit_errors(rdout_pos(idx),phenc_pos(idx),:)))/norm(time_series);
        
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),:)=[r2s_est;A;error;0;0;fit_numvals];
    else
        tTimeit=0;
        r2_fit_results(rdout_pos(idx),phenc_pos(idx),end)=fit_numvals;
    end
    
end

%% Function r2_ARL
    function [r2s_est,A,nonlin_part]=r2_ARL(time_series,M,alpha,te)
        
        temp1=time_series(M);
        
        S=alpha*(temp1(1,:)+4*temp1(2,:)+temp1(3,:));
        d=temp1(1,:)-temp1(3,:);
        
        p=sum(S.^2);
        q=sum(S.*d);
        r2s_est=(alpha*p + q)/(p + alpha*q);
        nonlin_part=exp(-te(:)*r2s_est);
        A=nonlin_part\time_series(:);
    end
end