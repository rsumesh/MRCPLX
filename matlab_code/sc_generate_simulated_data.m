%% Script to simulate T1, T2 and T2* data

% Author: Suryanarayana Umesh Rudrapatna
% e-mail: umeshrs at gmail.com
% Release date: Aug 29, 2015
% Reference: S Umesh Rudrapatna et. al., Improved estimation of MR relaxation parameters using complex-valued data
% MRCPLX  Copyright (C) {2015}
% This program comes with ABSOLUTELY NO WARRANTY
% Note: Toggle the flags for the if statements to generate desired data


% Common parameters for T1, T2 and T2* data generation

N=10;%Number of time series, 1000 used for data analyzed in the manuscript
A=1.0;%Proton density
op_dir='';% Directory where the simulated data is to be saved

%% T2* related
if 1
    snrs=[30 20 10];% SNRs for simulations
    noise_std=A./snrs;% Noise standard deviations
    ne_t2s=25;%Number of echoes
    te=2.5e-3;%TE gap for T2* acquisition
    T2s=[10 20 30 50]*1e-3;%Typical range of T2* values
    f_range=[-10 0 25 50 100];% Hz, frequency offsets for T2* simulations
    start_ph=linspace(0,2*pi,9);%Various starting phases
    start_ph=start_ph(1:end-1);
    
    tes=te*(1:ne_t2s);
    t2s_mag=exp(-tes(:)*transpose(1./T2s(:)));
    t2s_mag=repmat(t2s_mag,[1 1 N]);
    t2s_mag=permute(t2s_mag,[3 2 1]);
    szdata=size(t2s_mag);
    
    t2s_data=zeros([size(t2s_mag) length(f_range) length(start_ph) length(snrs)]);
    
    % Add phase and noise to the data
    for snr_idx=1:length(snrs)
        for f_idx=1:length(f_range)
            for ph_idx=1:length(start_ph);
                t2s_ph=start_ph(ph_idx) + 2*pi*f_range(f_idx)*tes(:);
                t2s_ph=reshape(t2s_ph,[1 1 length(t2s_ph)]);
                t2s_ph=repmat(t2s_ph,[szdata(1:2) 1]);
                t2s_cplx_data=t2s_mag.*exp(1i*t2s_ph);
                t2s_data(:,:,:,f_idx,ph_idx,snr_idx)=t2s_cplx_data+noise_std(snr_idx)*randn(szdata)+1i*noise_std(snr_idx)*randn(szdata);
            end
        end
    end
    save([op_dir 't2s_data_sim.mat'],'snrs','ne_t2s','te','T2s','f_range','start_ph','tes','t2s_data','noise_std');
end

%% T2 related
if 0
    snrs=[30 20 10];% SNR for simulations
    noise_std=A./snrs;% Noise standard deviations
    ne_t2=12;%Number of echoes
    te=12e-3;%TE for T2 acq
    T2=[20 40 60 80]*1e-3;%Typical range of T2 values
    f_range=0;% Hz, frequency offsets for T2* simulations
    start_ph=linspace(0,2*pi,9);%Various starting phases
    start_ph=start_ph(1:end-1);
    
    tes=te*(1:ne_t2);
    t2_mag=exp(-tes(:)*transpose(1./T2(:)));
    t2_mag=repmat(t2_mag,[1 1 N]);
    t2_mag=permute(t2_mag,[3 2 1]);
    szdata=size(t2_mag);
    
    t2_data=zeros([size(t2_mag) length(f_range) length(start_ph) length(snrs)]);
    
    % Add phase and noise to the data
    for snr_idx=1:length(snrs)
        for f_idx=1:length(f_range)
            for ph_idx=1:length(start_ph);
                t2_ph=start_ph(ph_idx) + 2*pi*f_range(f_idx)*tes(:);
                t2_ph=reshape(t2_ph,[1 1 length(t2_ph)]);
                t2_ph=repmat(t2_ph,[szdata(1:2) 1]);
                t2_cplx_data=t2_mag.*exp(1i*t2_ph);
                t2_data(:,:,:,f_idx,ph_idx,snr_idx)=t2_cplx_data+noise_std(snr_idx)*randn(szdata)+1i*noise_std(snr_idx)*randn(szdata);
            end
        end
    end
    save([op_dir 't2_data_sim.mat'],'snrs','ne_t2','te','T2','f_range','start_ph','tes','t2_data','noise_std');
end

%% T1 related
if 0
    snrs=[5 3 2];% SNR for simulations
    noise_std=A./snrs;% Noise standard deviations
    T1=[0.5 1.0 1.5 2.0 2.5];%Typical range of T1 values
    ti= 0.0560 + (0:148)*0.0480;%From our experimental settings
    bbya=2.0;% Parameter p in the article
    ne_t1=length(ti);% No. of data points
    start_ph=linspace(0,2*pi,9);%Various starting phases
    start_ph=start_ph(1:end-1);
    
    t1_mag=(1-(bbya*exp(-ti(:)*(bbya-1)*transpose(1./T1(:)))));
    t1_mag=repmat(t1_mag,[1 1 N]);
    t1_mag=permute(t1_mag,[3 2 1]);
    
    szdata=size(t1_mag);
    t1_data=zeros([szdata length(start_ph) length(snrs)]);
    
    % Add phase and noise to the data
    for snr_idx=1:length(snrs)
        for ph_idx=1:length(start_ph);
            t1_ph=start_ph(ph_idx);
            t1_cplx_data=t1_mag*exp(1i*t1_ph);
            t1_data(:,:,:,ph_idx,snr_idx)=t1_cplx_data+noise_std(snr_idx)*randn(szdata)+1i*noise_std(snr_idx)*randn(szdata);
        end
    end
    save([op_dir 't1_data_sim.mat'],'snrs','T1','ti','bbya','ne_t1','start_ph','t1_data','noise_std');
end