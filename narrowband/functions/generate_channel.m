function [G_cas, phi_all] = generate_channel(para)
%Generate cascaded channels
%  [G_cas, phi_all] = generate_channel(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   G_cas: cascaded channel for all users
%   phi_all: directions of the pathes between BS and STARS
%Date: 27/06/2023
%Author: Zhaolin Wang


% get the HITRAN matrix
HITRANparams = importdata('data_freq_abscoe.txt');
L = 4;

%% BS to STARS channel
% pathloss
D = 10;
path_loss_bs = getSpreadLoss(para.f, D) + getAbsLoss(para.f, D, HITRANparams );

G = 0;
phi_all = zeros(L,1);
for i = 1:L
    phi_all(i) = unifrnd(-pi/2, pi/2);
    alpha = sqrt(1/2)*(randn(1) + 1i*randn(1));
    b = steering_vector_ULA(phi_all(i), para.N);
    a = steering_vector_UPA(unifrnd(-pi/2, pi/2), unifrnd(-pi/2, pi/2), para.M_h, para.M_v);
    G = G + alpha*a*b';    
end





%% STAR-RIS to users channel
% pathloss
D = 3;

path_loss_su = getSpreadLoss(para.f, D) + getAbsLoss(para.f, D, HITRANparams );
path_loss = 10.^((-path_loss_bs - path_loss_su - para.noise_dB + para.Gt + para.Gr)/10);

h = zeros(para.M,para.K);
for k = 1:para.K
    h(:,k) = zeros(para.M,1);
    for i = 1:L
        alpha = sqrt(1/2)*(randn(1) + 1i*randn(1));
        h(:,k) = h(:,k) + alpha*steering_vector_UPA(unifrnd(-pi/2, pi/2), unifrnd(-pi/2, pi/2), para.M_h, para.M_v); 
    end
end


%% cascaded channel
G_cas = zeros(para.M, para.N, para.K);
for k = 1:para.K
    G_cas(:,:,k) = sqrt(path_loss)*diag(h(:,k)) * G; 
end

end
