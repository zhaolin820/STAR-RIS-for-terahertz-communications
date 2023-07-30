function [SE,R] = sum_rate(para, theta_t, theta_r, F_RF, T, F_BB, G)
%Calculate spectral efficiency achieved by hybrid beamforming
%  [SE,R] = sum_rate(para, theta_t, theta_r, F_RF, F_BB, G)
%Inputs:
%   para: structure of the initial parameters
%   theta_t, theta_r, F_RF, T, F_BB: obtainded beamformers
%   G: cascaded channel for all users
%Outputs:
%   SE: spectral efficiency
%   R: rate of each user
%Date: 27/07/2023
%Author: Zhaolin Wang

R = zeros(para.K,para.Mc);

for m = 1:para.Mc
    F_BB_m = F_BB(:,:,m);
    Tm = T(:,:,m);
    for k = 1:para.K
        
        if k <= para.K/2
            theta_i = theta_r;
        else
            theta_i = theta_t;
        end
        
        % calculate inter-user interference
        F_BB_inter_m = F_BB_m; F_BB_inter_m(:,k) = [];
        I_km = norm(theta_i.' * G(:,:,k, m) * F_RF * Tm * F_BB_inter_m)^2;
    
        % power of desired signal 
        p_k = abs(theta_i.' * G(:,:,k, m) * F_RF * Tm * F_BB_m(:,k))^2; 
    
        % rate for user k
        R(k,m) = log2( 1 + p_k / (I_km + 1) );
    end
end

% Spectral efficiency
SE = 1 / (para.Mc + 4) * sum(sum(R));

end

