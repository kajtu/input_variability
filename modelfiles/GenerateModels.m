% Compile model files into mexfiles
% This automatically creates the files: amai_avatar_corr.mex.. and simulate_avatar_corr.m

clear 
clear mex
 
addpath(genpath('..\Requirements' ))
run(['..' filesep 'Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab' filesep 'installAMICI.m'])

path = pwd;

amiwrap('avatar_belen','avatar_syms_belen',path); 

disp('Models generated')

