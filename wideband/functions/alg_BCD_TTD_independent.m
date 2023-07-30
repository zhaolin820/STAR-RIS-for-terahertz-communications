function [theta_t, theta_r, F, F_RF, T, F_BB, t, P, a, b] = alg_BCD_TTD_independent(para, G, rho, w, theta_t, theta_r, F_RF, T, F_BB, P, a, b, Psi, lambda, index_F_RF, index_T, t)
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
%   index_F_RF, index_T: indicator of the zero elements in F_RF and T
%   t: time delays of each TTD
%Outputs:
%   theta_t, theta_r, F_RF, T, F_BB: beamformers
%   t: time delays of each TTD
%   F, P, a, b: auxiliary variables
%Date: 27/07/2021
%Author: Zhaolin Wang

obj_pre = 100;
for i = 1:30

    % optimization of blocks
    [F, P, a, b, eta] = update_F(para, G, rho, w, theta_t, theta_r, F_RF, T, F_BB, P, a, b, Psi, lambda);
    [theta_t, theta_r] = update_theta(para, G, rho, F, P, theta_t, theta_r, lambda);
    [F_RF] = update_F_RF(para, rho, F, F_RF, T, F_BB, Psi, index_F_RF);
    [F_BB] = update_F_BB(para, rho, F, F_RF, T, Psi);
    [T, t] = update_T(para, rho, F, F_RF, F_BB, Psi, t, index_T);
    
    % calculate objective value
    penalty_1 = 0; penalty_2 = 0;
    for m = 1:para.Mc
        penalty_1 = penalty_1 + sum_square_abs( vec(F(:,:,m) - F_RF*T(:,:,m)*F_BB(:,:,m) + rho*Psi(:,:,m)) );
        for k = 1:para.K
            if k <= para.K/2
                theta_i = theta_r;
            else
                theta_i = theta_t;
            end
    
            penalty_2 = penalty_2 + sum_square_abs( P(k,:,m) - theta_i.'*G(:,:,k,m)*F(:,:,m) + rho*lambda(k,:,m) );
        end
    end

    obj = eta - 1/(2*rho) * (penalty_1 + penalty_2);
    
    % display
    [SE] = sum_rate(para, theta_t, theta_r, F_RF, T, F_BB, G);
    Pt = 0;
    for m = 1:para.Mc 
        Pt = Pt + norm(F_RF * T(:,:,m) * F_BB(:,:,m), 'fro')^2;
    end
    Pt_sum = para.Pc_TD_idp + 1/para.Mc*Pt + para.xi*SE; % power consumption
    EE = SE / Pt_sum;
    disp(['inner loop, i - ' num2str(i) ', eta - ' num2str(eta) ', penalty - ' num2str(penalty_1 + penalty_2)...
        ', SE - ' num2str(SE) ', EE - ' num2str(EE)]);
    
    % break or not?
    if abs(obj - obj_pre) / obj < 1e-3
        break
    end
    obj_pre = obj;
end

end

%% optimize auxiliary variables
function [F, P, a, b, eta] = update_F(para, G, rho, w, theta_t, theta_r, F_RF, T, F_BB, P, a, b, Psi, lambda)
a_t = a; b_t = b; P_t = P;
v = 1 / (para.Mc+4);
cvx_begin quiet
    % optimization variables
    variable F(para.N, para.K, para.Mc) complex
    variable P(para.K, para.K, para.Mc) complex
    variable r(para.K, para.Mc)
    variables eta a b

    % constraints
    eta <= 2*a_t/b_t * a - (a_t/b_t)^2 * b;
    Pt_sum = 0;
    for m = 1:para.Mc
        Pt_m = sum_square_abs(vec(F(:,:,m)));
        Pt_m <= para.Pt;
        Pt_sum = Pt_sum + Pt_m;
        for k = 1:para.K
            r(k,m) >= 0;
            
            P_km = P(k,:,m); P_km(k) = [];
            P_km_t = P_t(k,:,m); P_km_t(k) = [];
            I_k_m = sum_square_abs(P_km) + 1;
            I_k_m_t = sum_square_abs(P_km_t) + 1;
            
            gamma_mk = 2 * real( conj(P_t(k,k,m)) * P(k,k,m) ) / I_k_m_t - pow_abs(P_t(k,k,m)/I_k_m_t, 2) * I_k_m;
            power(2,r(k,m))-1 <= gamma_mk;
        end
    end
    a^2 <= v*sum(sum(r));
    w * ( 1/para.Mc*Pt_sum + para.xi * v*sum(sum(r)) ) + para.Pc_TD_idp <= b;

    % objective function
    penalty_1 = 0; penalty_2 = 0;
    for m = 1:para.Mc
        penalty_1 = penalty_1 + sum_square_abs( vec(F(:,:,m) - F_RF*T(:,:,m)*F_BB(:,:,m) + rho*Psi(:,:,m)) );
        for k = 1:para.K
            if k <= para.K/2
                theta_i = theta_r;
            else
                theta_i = theta_t;
            end
    
            penalty_2 = penalty_2 + sum_square_abs( P(k,:,m) - theta_i.'*G(:,:,k,m)*F(:,:,m) + rho*lambda(k,:,m) );
        end
    end

    obj = eta - 1/(2*rho) * (penalty_1 + penalty_2);

    maximize(obj);
cvx_end

end

%% optimiza STAR-RIS coefficients
function [theta_t, theta_r] = update_theta(para, G, rho, F, Xi, theta_t, theta_r, lambda)

Phi_t = zeros(para.M, para.M); Phi_r = zeros(para.M, para.M);
v_t = zeros(para.M, 1); v_r = zeros(para.M, 1);

for m = 1:para.Mc
    Fm = F(:,:,m);
    for k = 1:para.K
        G_mk = G(:,:,k,m);
        xi_mk = Xi(k,:,m);
        lambda_mk = lambda(k,:,m);
        if k <= para.K/2
            Phi_r = Phi_r + conj(G_mk)*conj(Fm)*Fm.'*G_mk.';
            v_r = v_r + conj(G_mk)*conj(Fm)*(xi_mk.' + rho*lambda_mk.');
        else
            Phi_t = Phi_t + conj(G_mk)*conj(Fm)*Fm.'*G_mk.';
            v_t = v_t + conj(G_mk)*conj(Fm)*(xi_mk.' + rho*lambda_mk.');
        end
    end
end

penalty_pre = 100;
for i = 1:50
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
    for m = 1:para.Mc
        for k = 1:para.K
            if k <= para.K/2
                theta_i = theta_r;
            else
                theta_i = theta_t;
            end
    
            penalty = penalty + sum_square_abs( Xi(k,:,m) - theta_i.'*G(:,:,k,m)*F(:,:,m) + rho*lambda(k,:,m) );
        end
    end
    
    % break or not?
    if abs(penalty - penalty_pre) / penalty < 1e-3
        break;
    end
    
    penalty_pre = penalty;  
end
end

%% optimize analog beamforming
function [F_RF] = update_F_RF(para, rho, F, F_RF, T, F_BB, Psi, index_F_RF)
A = 0; B = 0;
for m = 1:para.Mc
    A = A + T(:,:,m) * F_BB(:,:,m) * (T(:,:,m) * F_BB(:,:,m))';
    B = B + (F(:,:,m) + rho*Psi(:,:,m)) * (T(:,:,m) * F_BB(:,:,m))';
end

% BCD-based algorithm
penalty_pre = 100;
for n = 1:200
    for i = 1:para.N
        for j = 1:para.N_RF*para.N_T
            if index_F_RF(i,j) == 1
                Q = F_RF*A;
                q_ij = F_RF(i,j)*A(j,j) - Q(i,j) + B(i,j);
                F_RF(i,j) = q_ij / abs(q_ij);
            else
                F_RF(i,j) = 0;
            end
        end
    end
    
    % break or not?
    penalty = 0;
    for m = 1:para.Mc
        penalty = penalty + sum_square_abs( vec(F(:,:,m) - F_RF*T(:,:,m)*F_BB(:,:,m) + rho*Psi(:,:,m)) );
    end
    if abs(penalty - penalty_pre) / penalty < 1e-3
        break;
    end
    penalty_pre = penalty;
end
end

%% optimize digital beamforming
function [F_BB] = update_F_BB(para, rho, F, F_RF, T, Psi)
F_BB = zeros(para.N_RF, para.K, para.Mc);
for m = 1:para.Mc
    T_m = T(:,:,m);
    F_BB(:,:,m) = inv((F_RF*T_m)'*(F_RF*T_m)) * (F_RF*T_m)' * (F(:,:,m)+rho*Psi(:,:,m));
end

end

%% optimize TTD matrix
function [T, t] = update_T(para, rho, F, F_RF, F_BB, Psi, t, index_T)
C = F_RF'*F_RF;
D = zeros(para.K, para.K, para.Mc);
E = zeros(para.N_RF*para.N_T, para.N_RF, para.Mc);
for m = 1:para.Mc
    D(:,:,m) = F_BB(:,:,m)*F_BB(:,:,m)';
    E(:,:,m) = F_RF'*(F(:,:,m) + rho*Psi(:,:,m))*F_BB(:,:,m)';
end

% quasi-newton method
t_init = t;
options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
t = fminunc(@(t)(TTD_obj(t, para.fm_all, para.Mc, C, D, E, index_T)),t_init,options);

% calculate TTD matrix
t_diag = num2cell(t, 1);
t_diag = blkdiag(t_diag{:});

T = zeros(para.N_T * para.N_RF, para.N_RF, para.Mc);
for m = 1:para.Mc
    T_m = exp(1i*2*pi*para.fm_all(m)*t_diag); T_m(index_T ~= 1) = 0;
    T(:,:,m) = T_m;
end

end


function [f] = TTD_obj(t, fm_all, Mc, C, D, E, index_T)
t_diag = num2cell(t, 1);
t_diag = blkdiag(t_diag{:});

% Calculate objective f
f = 0;
for m = 1:Mc
    Tm = exp(1i*2*pi*fm_all(m)*t_diag); Tm(index_T ~= 1) = 0;
    f = f + real(trace(Tm'*C*Tm*D(:,:,m))) - 2*real(trace(Tm'*E(:,:,m)));
end

end

