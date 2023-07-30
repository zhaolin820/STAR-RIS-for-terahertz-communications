function [theta_t, theta_r, F, F_RF, F_BB, P, a, b] = alg_BCD_independent(para, G, rho, w, theta_t, theta_r, F_RF, F_BB, P, a, b, Psi, lambda)
%BCD-based algorithm: inner loop 
%  alg_BCD_independent(para, G, rho, w, theta_t, theta_r, F_RF, F_BB, P, a, b, Psi, lambda)
%Inputs:
%   para: structure of the initial parameters
%   G: cascaded channel for all users
%   rho: panalty factor
%   phi_all: directions of the pathes between BS and STARS
%   w: w=1 for maximizing EE; w=0 for maximizing SE
%   theta_t, theta_r, F_RF, F_BB: beamformers
%   P, a, b: auxiliary variables
%   Psi, lambda: dual variables
%Outputs:
%   theta_t, theta_r, F_RF, F_BB: beamformers
%   F, P, a, b: auxiliary variables
%Date: 27/06/2021
%Author: Zhaolin Wang

obj_pre = 100;
for i = 1:30

    % optimization of blocks
    [F, P, a, b, eta] = update_F(para, G, rho, w, theta_t, theta_r, F_RF, F_BB, P, a, b, Psi, lambda);
    [theta_t, theta_r] = update_theta(para, G, rho, F, P, theta_t, theta_r, lambda);
    [F_RF] = update_F_RF(para, rho, F, F_RF, F_BB, Psi);
    [F_BB] = update_F_BB(rho, F, F_RF, Psi);
    
    % calculate objective value
    penalty_1 = sum_square_abs( vec(F - F_RF*F_BB + rho*Psi) );
    penalty_2 = 0;
    for k = 1:para.K
        if k <= para.K/2
            theta_i = theta_r;
        else
            theta_i = theta_t;
        end

        penalty_2 = penalty_2 + sum_square_abs( P(k,:) - theta_i.'*G(:,:,k)*F + rho*lambda(k,:) );
    end
    obj = eta - 1/(2*rho) * (penalty_1 + penalty_2);

    % display
    [SE] = sum_rate(para, theta_t, theta_r, F_RF, F_BB, G);
    EE = SE / (para.Pc_HB_idp + norm(F_RF*F_BB, 'fro')^2 + para.xi*SE);
    disp(['inner loop, i - ' num2str(i) ', eta - ' num2str(eta) ', penalty - ' num2str(penalty_1+penalty_2)...
        ', SE - ' num2str(SE) ', EE - ' num2str(EE)]);
    
    % break or not?
    if abs(obj - obj_pre) / obj < 1e-3
        break
    end
    obj_pre = obj;
end

end

%% optimize auxiliary variables
function [F, P, a, b, eta] = update_F(para, G, rho, w, theta_t, theta_r, F_RF, F_BB, P, a, b, Psi, lambda)
a_t = a; b_t = b; P_t = P;

cvx_begin quiet
    % optimization variables
    variable F(para.N, para.K) complex
    variable P(para.K, para.K) complex
    variable r(para.K, 1)
    variables eta a b

    % constraints
    eta <= 2*a_t/b_t * a - (a_t/b_t)^2 * b;
    for k = 1:para.K
        r(k) >= 0;

        I_k = 0;
        I_k_t = 0;
        for j = 1:para.K
            if j ~= k
                I_k = I_k + pow_abs(P(k,j),2);
                I_k_t = I_k_t + pow_abs(P_t(k,j),2);
            end
        end
        I_k = I_k + 1;
        I_k_t = I_k_t + 1;
        gamma_k = 2 * real( conj(P_t(k,k)) * P(k,k) ) / I_k_t - pow_abs(P_t(k,k)/I_k_t, 2) * I_k;
        power(2,r(k))-1 <= gamma_k;
    end
    a^2 <= sum(r);
    w * ( sum_square_abs(vec(F)) + para.xi * sum(r) ) + para.Pc_HB_idp <= b;
    sum_square_abs(vec(F)) <= para.Pt;

    % objective function
    penalty_1 = sum_square_abs( vec(F - F_RF*F_BB + rho*Psi) );
    penalty_2 = 0;
    for k = 1:para.K
        if k <= para.K/2
            theta_i = theta_r;
        else
            theta_i = theta_t;
        end
        penalty_2 = penalty_2 + sum_square_abs( P(k,:) - theta_i.'*G(:,:,k)*F + rho*lambda(k,:) );
    end

    obj = eta - 1/(2*rho) * (penalty_1 + penalty_2);

    maximize(obj);
cvx_end

end

%% optimize STAR-RIS coefficients
function [theta_t, theta_r] = update_theta(para, G, rho, F, P, theta_t, theta_r, lambda)

Phi_t = zeros(para.M, para.M); Phi_r = zeros(para.M, para.M);
v_t = zeros(para.M, 1); v_r = zeros(para.M, 1);

for k = 1:para.K
    G_k = G(:,:,k);
    p_k = P(k,:);
    lambda_k = lambda(k,:);
    if k <= para.K/2      
        Phi_r = Phi_r + conj(G_k)*conj(F)*F.'*G_k.';
        v_r = v_r + conj(G_k)*conj(F)*(p_k.' + rho*lambda_k.');
    else
        Phi_t = Phi_t + conj(G_k)*conj(F)*F.'*G_k.';
        v_t = v_t + conj(G_k)*conj(F)*(p_k.' + rho*lambda_k.');
    end
end

penalty_pre = 100;
while 1
    for m = 1:para.M
        c_t = Phi_t(m,m); c_r = Phi_r(m,m);
        t = Phi_t*theta_t; r = Phi_r*theta_r;
        d_t = Phi_t(m,m) * theta_t(m) - t(m) + v_t(m);
        d_r = Phi_r(m,m) * theta_r(m) - r(m) + v_r(m);

        % calculate optimal phase
        phi_t = angle(d_t);
        phi_r = angle(d_r);

        % calculate optimal amplitude
        fun = @(x)( c_t*(sin(x))^2 + c_r*(cos(x))^2 - 2*abs(d_t)*sin(x) - 2*abs(d_r)*cos(x));
        x = fminbnd(fun,0,pi/2);
        beta_t = sin(x); 
        beta_r = cos(x);     

        theta_t(m) = beta_t*exp(1i*phi_t);
        theta_r(m) = beta_r*exp(1i*phi_r);
    end
    
    % calculate penalty value
    penalty = 0;
    for k = 1:para.K
        if k <= para.K/2
            theta_i = theta_r;
        else
            theta_i = theta_t;
        end

        penalty = penalty + sum_square_abs( P(k,:) - theta_i.'*G(:,:,k)*F + rho*lambda(k,:));
    end
    
    % break or not?
    if abs(penalty - penalty_pre) / penalty < 1e-3
        break;
    end
    
    penalty_pre = penalty;  
end
end

%% optimize analog beamforming
function [F_RF] = update_F_RF(para, rho, F, F_RF, F_BB, Psi)
A = F_BB*F_BB'; B = (F + rho*Psi)*F_BB';
penalty_pre = 100;
for n = 1:200
    for i = 1:para.N
        for j = 1:para.N_RF
            Q = F_RF*A;
            q_ij = F_RF(i,j)*A(j,j) - Q(i,j) + B(i,j);
            F_RF(i,j) = q_ij / abs(q_ij);
        end
    end
    
    % break or not?
    penalty = norm(F - F_RF*F_BB + rho*Psi, 'fro');
    if abs(penalty - penalty_pre) / penalty < 1e-3
        break;
    end
    penalty_pre = penalty;
end
end

%% optimize digital beamforming
function [F_BB] = update_F_BB(rho, F, F_RF, Psi)

F_BB = inv(F_RF'*F_RF) * F_RF' * (F+rho*Psi);

end