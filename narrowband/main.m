clc
clear all
close all

addpath("functions\");
cvx_solver mosek

%% parameters and channels
para = para_init();
[G, phi_all] = generate_channel(para);

%% optimization
w = 1; % w=1 for maximizing EE; w=0 for maximizing SE;
[theta_t, theta_r, F_RF, F_BB] = alg_PDD_independent(para, G, phi_all, w);
[SE,R] = sum_rate(para, theta_t, theta_r, F_RF, F_BB, G);
EE = SE / (para.Pc_HB_idp + norm(F_RF*F_BB, 'fro')^2 + para.xi*SE);
