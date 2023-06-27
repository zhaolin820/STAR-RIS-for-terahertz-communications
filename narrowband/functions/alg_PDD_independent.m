function [theta_t, theta_r, F_RF, F_BB] = alg_PDD_independent(para, G, phi_all, w)
%PDD-based algorithm: outer loop 
%  [theta_t, theta_r, F_RF, F_BB] = alg_PDD_independent(para, G, phi_all, w)
%Inputs:
%   para: structure of the initial parameters
%   G: cascaded channel for all users
%   phi_all: directions of the pathes between BS and STARS
%   w: w=1 for maximizing EE; w=0 for maximizing SE
%Outputs:
%   theta_t, theta_r, F_RF, F_BB: obtainded beamformers
%Date: 27/06/2021
%Author: Zhaolin Wang

rho = 1e3; % penalty term
c = 0.6; % reduction factor of rho
epsilon = 1e-3; % convergence criteria

%% initialization

% initialize beamformers and STARS coefficients
theta_t = randn(para.M, 1) + 1i * randn(para.M, 1); theta_t = sqrt(0.5) * theta_t ./ abs(theta_t);
theta_r = randn(para.M, 1) + 1i * randn(para.M, 1); theta_r = sqrt(0.5) * theta_r ./ abs(theta_r);

F = randn(para.N, para.K) + 1i * randn(para.N, para.N_RF); 
F = sqrt(para.Pt)* F ./ norm(F, 'fro');

F_RF = zeros(para.N, 4);
for i = 1:4
    F_RF(:,i) = steering_vector_ULA(phi_all(i), para.N);
end


F_BB = inv(F_RF'*F_RF) * F_RF' * F;

[SE,~] = sum_rate_full_digital(para, theta_t, theta_r, F, G);

% initialize auxiliary variables
a = sqrt(SE);
b = w * (norm(F, 'fro')^2 + para.xi*SE) + para.Pc_HB_idp; % power consumption

P = zeros(para.K, para.K);
for k = 1:para.K
    if k <= para.K/2
        theta_i = theta_r;
    else
        theta_i = theta_t;
    end
    p_k =  theta_i.' * G(:,:,k) * F_RF * F_BB;
    P(k,:) = p_k;
end

% initialize dual variables
Psi = zeros(para.N, para.K);
lambda = zeros(para.K, para.K);


%% PDD algorthm
eta = 10;
for i = 1:40

    % optimize AL problem
    [theta_t, theta_r, F, F_RF, F_BB, P, a, b] = alg_BCD_independent(para, G, rho, w, theta_t, theta_r, F_RF, F_BB, P, a, b, Psi, lambda);
    
    [h] = constraint_violation(para, G, theta_t, theta_r, F, F_RF, F_BB, P); 
    disp(['%%%%%%%%%%%%%%%%%%%%%% Outter loop i - ' num2str(i) ', h - ' num2str(h) ' %%%%%%%%%%%%%%%%%%%%%']);

    
    if h < epsilon % algorithm converged
        break; 
    end

    if h <= eta    
        % update dual variables
        Psi = Psi + 1/rho * ( F - F_RF*F_BB );
        for k = 1:para.K
            if k <= para.K/2
                theta_i = theta_r;
            else
                theta_i = theta_t;
            end
            lambda(k,:) = lambda(k,:) + 1/rho * ( P(k,:) -  theta_i.'*G(:,:,k)*F);
        end
    else
        % update penalty term
        rho = c*rho; 
        disp(['rho - ' num2str(rho)]);
    end

    eta = 0.7*h;

end

end


function [h] = constraint_violation(para, G, theta_t, theta_r, F, F_RF, F_BB, Xi)
    a = norm(F - F_RF*F_BB, "inf");
    
    buffer = zeros(para.K, 1);
    for k = 1:para.K
        if k <= para.K/2
            theta_i = theta_r;
        else
            theta_i = theta_t;
        end
        buffer(k) = max( abs( Xi(k,:) -  theta_i.' * G(:,:,k) * F) );
    end
    b = max(buffer);
    
    h = max([a,b]);
end