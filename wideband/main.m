clc
clear all
close all

addpath("functions\");
cvx_solver mosek

rng(1);
%% parameters and channels
para = para_init();
[G, phi_all] = generate_channel(para);


%% optimization
w = 1; % w=1 for maximizing EE; w=0 for maximizing SE;
[theta_t, theta_r, F_RF, T, F_BB] = alg_PDD_TTD_independent(para, G, phi_all, w);
[SE, R] = sum_rate(para, theta_t, theta_r, F_RF, T, F_BB, G);
% Energy efficiency
Pt = 0;
for m = 1:para.Mc
    Pt = Pt + norm(F_RF * T(:,:,m) * F_BB(:,:,m), 'fro')^2;
end
P = para.Pc_TD_idp + 1/para.Mc*Pt + para.xi*SE; % power consumption
EE = SE / P;