function [para] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [values] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 27/07/2021
%Author: Zhaolin Wang

para.fc = 1e11; % carrier frequency 300 GHz
B = 1e10; % 10 GHz system bandwidth 

para.Mc = 10; % number of subcarriers

m = 1:para.Mc;
para.fm_all =  para.fc + B*(2*m-1-para.Mc) / (2*para.Mc); % subcarrier frequency

para.noise_dB = -174; % noise power in dBm/Hz
para.noise_dB = para.noise_dB + 10*log10(B/10) - 30; % noise power in dBW

para.Gt = 25; % dBi, transmit antenna gain
para.Gr = 25; % dBi, receive antenna gain

para.N = 128; % number of antennas
para.N_RF = 4; % number of RF chains
para.N_T = 8; % number of TTDs connected to each RF chain

para.M_h = 6;
para.M_v = 6;
para.M = para.M_h*para.M_v; % reflecting elements at RIS

para.Pt = 10^(20/10); % dBW overall transmit power

para.L = 4; % number of paths between BS and STARS
para.Lk = 4; % number of paths between STARS and each user
para.K = 4; % user number


P_UE = 100; % 100 mW, user power 
P_STAR_independent = para.M * 1/2*(7+2*8)*0.33 + 1e4; % mW, STAR-RIS power
P_STAR_coupled = para.M * 1/2*(7+8+1)*0.33 + 1e4; % mW, STAR-RIS power
P_RIS = para.M * 1/2*(8)*0.33 + 1e4;
P_BB = 300; % 300 mW, baseband processing power 
P_RF = 200; % 200 mW, RF chain power
P_PS = 20; % 20 mW, PS power
P_BS = 3e3; % 3 W, BS power
P_TTD = 100; % 100 mW, TTD power

para.Pc_TD_idp =1e-3 * (P_BS + P_BB + para.N_RF * P_RF + para.N_RF * para.N_T * P_TTD + para.N_RF * para.N * P_PS + P_STAR_independent + para.K*P_UE);
para.Pc_TD_coup =1e-3 * (P_BS + P_BB + para.N_RF * P_RF + para.N_RF * para.N_T * P_TTD + para.N_RF * para.N * P_PS + P_STAR_coupled + para.K*P_UE); 

para.xi = 1e-3 * 100; % 10 mW/(bit/s/Hz), rate-dependent power consumption factor.

end

