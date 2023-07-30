function [SE,R] = sum_rate_full_digital(para, theta_t, theta_r, F, G)
%Calculate spectral efficiency achieved by fully-digital beamforming
%  [SE,R] = sum_rate(para, theta_t, theta_r, F_RF, F_BB, G)
%Inputs:
%   para: structure of the initial parameters
%   theta_t, theta_r, F: obtainded beamformers
%   G: cascaded channel for all users
%   theta_v: elevation angle
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
    F_inter = F; F_inter(:,k) = [];
    I_k = norm(theta_i.' * G(:,:,k) * F_inter)^2;

    % power of desired signal 
    p_k = abs(theta_i.' * G(:,:,k) * F(:,k))^2; 

    % rate for user k
    R(k) = log2( 1 + p_k / (I_k + 1) );
end

% Spectral efficiency
SE = sum(R);

end

