function [theta_t, theta_r, F_RF, T, F_BB] = alg_PDD_TTD_independent(para, G, phi_all, w)
%PDD-based algorithm: outer loop 
%  [theta_t, theta_r, F_RF, F_BB] = alg_PDD_independent(para, G, phi_all, w)
%Inputs:
%   para: structure of the initial parameters
%   G: cascaded channel for all users
%   phi_all: directions of the pathes between BS and STARS
%   w: w=1 for maximizing EE; w=0 for maximizing SE
%Outputs:
%   theta_t, theta_r, F_RF, T, F_BB: obtainded beamformers
%Date: 27/07/2021
%Author: Zhaolin Wang

rho = 1e3; % penalty term
c = 0.6; % reduction factor of rho
epsilon = 5e-2; % convergence criteria

%% initialization

% initialize beamformers and STARS coefficients
theta_t = randn(para.M, 1) + 1i * randn(para.M, 1); theta_t = sqrt(0.5) * theta_t ./ abs(theta_t);
theta_r = randn(para.M, 1) + 1i * randn(para.M, 1); theta_r = sqrt(0.5) * theta_r ./ abs(theta_r);

F_RF = zeros(para.N, para.N_RF*para.N_T);
t = zeros(para.N_T, para.N_RF);
index_T = zeros(para.N_T * para.N_RF, para.N_RF);
N_sub = para.N / para.N_T;
for l = 1:4
    theta_l = sin(phi_all(l));
    f = steering_vector_ULA(phi_all(l), para.N, para.fc, para.fc);
    F_RF_l = zeros(para.N, para.N_T);
    t_l = zeros(para.N_T, 1);
    
    for i = 1:para.N_T
        F_RF_l((i-1)*N_sub+1:i*N_sub, i) = f((i-1)*N_sub+1:i*N_sub) * exp(-1i*pi*(i-1)*N_sub*theta_l);
    
        if theta_l >= 0
            t_l(i) = (i-1) * (N_sub*theta_l/2) * 1/para.fc;
        else
            t_l(i) = (para.N_T-1) * abs(N_sub*theta_l/2) * 1/para.fc + (i-1)* (N_sub*theta_l/2) * 1/para.fc;
        end
    end
    
    t(:, l) = t_l;   
    F_RF(:, (l-1)*para.N_T+1: l*para.N_T) = F_RF_l;
    index_T((l-1)*para.N_T+1: l*para.N_T, l) = ones(para.N_T, 1);
end
index_F_RF = F_RF ~= 0;

t_diag = num2cell(t, 1);
t_diag = blkdiag(t_diag{:});

T = zeros(para.N_T * para.N_RF, para.N_RF, para.Mc);
F_BB = zeros(para.N_RF, para.K, para.Mc);
F = zeros(para.N, para.K, para.Mc);
for m = 1:para.Mc
    T_m = exp(1i*2*pi*para.fm_all(m)*t_diag); T_m(index_T ~= 1) = 0;
    T(:,:,m) = T_m;

    F(:,:,m) = randn(para.N, para.K) + 1i * randn(para.N, para.K); 
    F(:,:,m) = sqrt(para.Pt)* F(:,:,m) ./ norm(F(:,:,m), 'fro');

    F_BB(:,:,m) = inv((F_RF*T_m)'*(F_RF*T_m)) * (F_RF*T_m)' * F(:,:,m);
end


% initialize auxiliary variables
[SE,~] = sum_rate_full_digital(para, theta_t, theta_r, F, G);
a = sqrt(SE);
Pt_used = 0;
for m = 1:para.Mc
    Pt_used = Pt_used + norm(F_RF * T(:,:,m) * F_BB(:,:,m), 'fro')^2;
end
b = w * (1/para.Mc*Pt_used + para.xi*SE) + para.Pc_TD_idp; % power consumption


P = zeros(para.K, para.K, para.Mc);
for m = 1:para.Mc
    for k = 1:para.K
        if k <= para.K/2
            theta_i = theta_r;
        else
            theta_i = theta_t;
        end
        p_km =  theta_i.' * G(:,:,k,m) * F_RF * T(:,:,m) * F_BB(:,k,m);
        P(k,:,m) = p_km;
    end
end

% dual variables
Psi = zeros(para.N, para.K, para.Mc);
lambda = zeros(para.K, para.K, para.Mc);


%% PDD algorthm
varepsilon = 10;
for i = 1:40

    % optimize AL problem
    [theta_t, theta_r, F, F_RF, T, F_BB, t, P, a, b] = alg_BCD_TTD_independent(para, G, rho, w, theta_t, theta_r, F_RF, T, F_BB, P, a, b, Psi, lambda, index_F_RF, index_T, t);
    
    [h] = constraint_violation(para, G, theta_t, theta_r, F, F_RF, T, F_BB, P); 
    disp(['%%%%%%%%%%%%%%%%%%%%%% Outter loop i - ' num2str(i) ', h - ' num2str(h) ' %%%%%%%%%%%%%%%%%%%%%']);
    
    if h < epsilon % algorithm converged
        break; 
    end

    if h <= varepsilon    
        % update dual variables
        for m = 1:para.Mc
            Psi(:,:,m) = Psi(:,:,m) + 1/rho * ( F(:,:,m) - F_RF*T(:,:,m)*F_BB(:,:,m) );
            for k = 1:para.K
                if k <= para.K/2
                    theta_i = theta_r;
                else
                    theta_i = theta_t;
                end
                lambda(k,:,m) = lambda(k,:,m) + 1/rho * ( P(k,:,m) -  theta_i.'*G(:,:,k,m)*F(:,:,m));
            end
        end
        disp('update dual variables');
    else
        % update penalty term
        rho = c*rho; 
        disp(['update penalty factor, rho - ' num2str(rho)]);
    end

    varepsilon = 0.7*h;

end

end


function [h] = constraint_violation(para, G, theta_t, theta_r, F, F_RF, T, F_BB, Xi)
    buffer_1 = zeros(para.Mc, 1);
    buffer_2 = zeros(para.Mc, para.K);
    for m = 1:para.Mc
        buffer_1(m) = norm( F(:,:,m) - F_RF*T(:,:,m)*F_BB(:,:,m),"inf" );
        for k = 1:para.K
            if k <= para.K/2
                theta_i = theta_r;
            else
                theta_i = theta_t;
            end
    
            buffer_2(m,k) = norm( Xi(k,:,m) - theta_i.'*G(:,:,k,m)*F(:,:,m),"inf" );
        end
    end

    a = max(buffer_1); b = max(max(buffer_2));
    h = max([a,b]);
end