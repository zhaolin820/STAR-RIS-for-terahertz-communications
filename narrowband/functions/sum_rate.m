function [SE,R] = sum_rate(para, theta_t, theta_r, F_RF, F_BB, G)
%Calculate spectral efficiency achieved by hybrid beamforming
%  [SE,R] = sum_rate(para, theta_t, theta_r, F_RF, F_BB, G)
%Inputs:
%   para: structure of the initial parameters
%   theta_t, theta_r, F_RF, F_BB: obtainded beamformers
%   G: cascaded channel for all users
%Outputs:
%   SE: spectral efficiency
%   R: rate of each user
%Date: 27/06/2023
%Author: Zhaolin Wang

R = zeros(para.K,1);

for k = 1:para.K
    
    if k <= para.K/2
        theta_i = theta_r;
    else
        theta_i = theta_t;
    end

    % calculate inter-user interference
    F_bb_inter = F_BB; F_bb_inter(:,k) = [];
    I_k = norm(theta_i.' * G(:,:,k) * F_RF * F_bb_inter)^2;

    % power of desired signal 
    p_k = abs(theta_i.' * G(:,:,k) * F_RF * F_BB(:,k))^2; 

    % rate for user k
    R(k) = log2( 1 + p_k / (I_k + 1) );
end

% Spectral efficiency
SE = sum(R);


end

