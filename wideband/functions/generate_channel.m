function [G_cas, phi_all] = generate_channel(para)
%Generate cascaded channels
%  [G_cas, phi_all] = generate_channel(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   G_cas: cascaded channel for all users
%   phi_all: directions of the pathes between BS and STARS
%Date: 27/07/2023
%Author: Zhaolin Wang

HITRANparams = importdata('data_freq_abscoe.txt');

%% BS to STARS channel

phi_all = rand(para.L,1)*pi-pi/2;
phi_2_all = rand(para.L,1)*pi-pi/2;
phi_3_all = rand(para.L,1)*pi-pi/2;
G = zeros(para.M, para.N, para.Mc);
for m = 1:para.Mc
    fm = para.fm_all(m);
    for i = 1:para.L
            alpha = sqrt(1/2)*(randn(1) + 1i*randn(1));
            b = steering_vector_ULA(phi_all(i), para.N, fm, para.fc);
            a = steering_vector_UPA(phi_2_all(i), phi_3_all(i), para.M_h, para.M_v, fm, para.fc);
            G(:,:,m) = G(:,:,m) + alpha*a*b';
    end
end

%% STARS to user channel
phi_4_all = rand(para.Lk,para.K)*pi-pi/2;
phi_5_all = rand(para.Lk,para.K)*pi-pi/2;

h = zeros(para.M,para.K, para.Mc);
for m = 1:para.Mc
    fm = para.fm_all(m);
    for k = 1:para.K
        for i = 1:para.Lk
            alpha = sqrt(1/2)*(randn(1) + 1i*randn(1));
            h(:,k,m) = h(:,k,m) + alpha*steering_vector_UPA(phi_4_all(i, k), phi_5_all(i, k), para.M_h, para.M_v, fm, para.fc);
        end
    end
end

%% cascaded channel
G_cas = zeros(para.M, para.N, para.K, para.Mc);
for m = 1:para.Mc
    fm = para.fm_all(m);
    for k = 1:para.K
        % frequency-dependent pathloss
        D = 10; path_loss_bs = getSpreadLoss(fm, D) + getAbsLoss(fm, D, HITRANparams );
        D = 3; path_loss_su = getSpreadLoss(fm, D) + getAbsLoss(fm, D, HITRANparams );
        path_loss = 10.^((-path_loss_bs - path_loss_su - para.noise_dB + para.Gt + para.Gr)/10);

        G_cas(:,:,k,m) = sqrt(path_loss)*diag(h(:,k,m)) * G(:,:,m);
    end
end




end