% Performing sensitivity analysis of the cardiovascular model by changing
% the values of one parameter at a time and observing the resulting model
% simulations.

close all
clear

%% Setup amici toolbox and make sure the simulation and model scripts are accessible
addpath(genpath('.'))
run(['Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab' filesep 'installAMICI.m'])

%% Set parameter values and parameter names
% select the parameter values for one of the subjects as the basis of the sensitivity analysis
parameterValues = [6.49538991133425	0.123134510159502	3.07934797941818	0.226561358768333	1.49082946131413	0.135811466173948	0.0816420562775852	0.000100017413547206	0.000345919441150104	0.000789742456908828	11.6970990247282	0.0799998566418412	0.00283178903696888	0.182079062245928	0.427059946089428	0.110112595957737	0.352054081810405	0.669571952186817	1.27539069941550	20.0629196863837	33.2934275929550	0.839043306773090	-0.0402335002653342]; 
paramNames = {'A_ao'	'Caa'	'EOA_av'	'Emax_LA'	'Emax_LV'	'Emin_LA'	'Emin_LV'	'Lao'	'Lav'	'Lmv'	'Ppu'	'Rao'	'Rmv'	'k_diast_LA'	'k_diast_LV'	'k_syst_LA'	'k_syst_LV'	'm1_LA'	'm1_LV'	'm2_LA'	'm2_LV'	'onset_LA'	'onset_LV'};

minparamvalues=[4.41849379568160
0.0275485949608387
2.13288789167579
0.0852027846707382
1.33304046890193
0.0646241380314971
0.0400729917213614
0.000100000000026120
8.00000000009790e-05
0.000285961700059100
6.95841430476510
0.0578739295668872
0.00188884462657160
0.0932955945551530
0.360598641989237
0.0550024247019690
0.193212594685934
0.660252266241282
1.12345586545891
6.55121826243352
21.3796151149423
0.527575265448606
-0.110314979859732];
maxparamvalues=[10.5863592420565
0.200160482566724
4.44641767001152
0.296631011415156
3.13942105646184
0.159970260987769
0.152059042184548
0.000851401173319472
0.000593348214247623
0.000999750065213153
11.9926173013935
0.0799999995821364
0.00750177702767131
0.359245782029224
0.519281271005727
0.218379480789827
0.536951157751445
2.62932361009881
2.18630800388697
26.1665539143252
41.6669395590998
0.996782656869929
0.00704945017940273];

% Set values for the constant paramters
constantsNames = {'Cpvc' 'Rpu' 'Rpv' 'Lpv' 'Rtot' 'Ctot' 'Ks_LA' 'Ks_LV' 'V0_LA' 'V0_LV' 'RLAvisc' 'RLVvisc' 'Raa' 'Rpc' 'Rpvc' 'T' 'rho_blood' 'norm_factor_LA' 'norm_factor_LV'};
rho_blood=1.06; 
Rpu=0.01;
Rpvc=0.01;
Cpvc=4;
Lpv=0.0005;
Rpv=0.002;
Ks_LA=10e-9;  
Ks_LV=4e-9;   
V0_LA=3;    
V0_LV=10;         
RLAvisc=0.0001;
RLVvisc=0.0001; 
Raa=0.01;
Rpc=0.01;

qAV = [21.446836	109.42274	388.55142	500.02371	506.11694	486.63501	466.08621	433.94769	400.68817	353.18576	306.03964	249.86032	204.59236	144.91663	72.496414	12.247868	-4.9297948	-2.6157029	9.7481976	26.605375	38.033092	40.376442	33.297287	27.178595	12.398558	5.7927256	-4.9340363	-9.7044554	-8.8095455	-3.1921296	-1.8886228	-1.3096753	-0.22389048	4.2830038	11.544492	23.199047	20.106979	9.6216478	-5.9584904	-5.1960468];
time = [0	0.0240000000000000	0.0480000000000000	0.0730000000000000	0.0970000000000000	0.121000000000000	0.145000000000000	0.169000000000000	0.193000000000000	0.218000000000000	0.242000000000000	0.266000000000000	0.290000000000000	0.314000000000000	0.339000000000000	0.363000000000000	0.387000000000000	0.411000000000000	0.435000000000000	0.459000000000000	0.484000000000000	0.508000000000000	0.532000000000000	0.556000000000000	0.580000000000000	0.605000000000000	0.629000000000000	0.653000000000000	0.677000000000000	0.701000000000000	0.725000000000000	0.750000000000000	0.774000000000000	0.798000000000000	0.822000000000000	0.846000000000000	0.871000000000000	0.895000000000000	0.919000000000000	0.943000000000000];
T = time(end) + mean(diff(time));
HR = 60/T;

% calculate parameters
SBP=120;
DBP=70;
PP=SBP-DBP;         %Peak pressure
MAP=DBP+(1/3+HR*0.0012)*(SBP-DBP); % Mean arterial pressure, takes HR in consideration
Ctot=mean(qAV)/PP;  %Total compliance of the system
Rtot=MAP/mean(qAV);     %Total resistance of the system

%% Set simulation settings
step = 0.001; % simulation time step
simtime = 0:step:T; % time vector to be simulated (for each heartbeat)
options.x0 = [7.322,7.751,94.333,7.751,14.789,7.751,7.275,7.751,6.965,1,1]'; % initial conditions for the simulation

% set indexes needed in the simulation function simulate_avatarHEALTH
indexes.T = strcmp('T',constantsNames);
indexes.k_syst_LV = strcmp('k_syst_LV',paramNames);
indexes.k_syst_LA = strcmp('k_syst_LA',paramNames);
indexes.k_diast_LA = strcmp('k_diast_LA',paramNames);
indexes.k_diast_LV = strcmp('k_diast_LV',paramNames);
indexes.m1_LA = strcmp('m1_LA',paramNames);
indexes.m1_LV = strcmp('m1_LV',paramNames);
indexes.m2_LA = strcmp('m2_LA',paramNames);
indexes.m2_LV = strcmp('m2_LV',paramNames);
indexes.onset_LV = strcmp('onset_LV',paramNames);
indexes.onset_LA = strcmp('onset_LA',paramNames);

% make sure onset is ok
if parameterValues(indexes.onset_LV) <0
    parameterValues(indexes.onset_LV) =parameterValues(indexes.onset_LV) + T;
end
if parameterValues(indexes.onset_LA) <0
    parameterValues(indexes.onset_LA) =parameterValues(indexes.onset_LA) + T;
end

% Calculate normalizing factors for elastance function based on the parameters
norm_factor_LV = calc_norm_factor(T,parameterValues(indexes.k_syst_LV),parameterValues(indexes.k_diast_LV),parameterValues(indexes.m1_LV),parameterValues(indexes.m2_LV));
norm_factor_LA = calc_norm_factor(T,parameterValues(indexes.k_syst_LA),parameterValues(indexes.k_diast_LA),parameterValues(indexes.m1_LA),parameterValues(indexes.m2_LA));

constants = [Cpvc Rpu Rpv Lpv Rtot Ctot Ks_LA Ks_LV V0_LA V0_LV RLAvisc RLVvisc Raa Rpc Rpvc T rho_blood norm_factor_LA norm_factor_LV];
constants = double(constants);
parameterValues = double(parameterValues);
parameterValuesOriginal = parameterValues;

%% Simulate the model
disp('Simulating...')
% simulate_avatarbelen creates repeated simulations of one heartbeat and
% sends out the last simulated heartbeat (simLast) (assumed to be at steady state).
% The heartbeat simulations are created with the basic simulation function
% simulate_avatar_belen that is generated by the AMICI toolbox.

numsims = 20; % number of parameter values to simualte for each investigated parameter
numHeartBeats = 20; %number of heartbeats to simulate until steady state is reached
extrarange = 0.20; % how far away from the min and max values the two more extreme simulations should be

%% simulate range of EOA values
eoarange = linspace(minparamvalues(strcmp(paramNames,'EOA_av')),maxparamvalues(strcmp(paramNames,'EOA_av')) ,numsims);
eoarange = [minparamvalues(strcmp(paramNames,'EOA_av')) - extrarange*minparamvalues(strcmp(paramNames,'EOA_av')),...
    eoarange,...
    maxparamvalues(strcmp(paramNames,'EOA_av')) + extrarange*maxparamvalues(strcmp(paramNames,'EOA_av'))];
warning 'off'
pvals = parameterValuesOriginal;
eoasims = cell(size(eoarange));
for i = 1:length(eoarange)
    pvals(strcmp(paramNames,'EOA_av')) = eoarange(i);
    [~,eoasims{i}] = simulate_avatarbelen(pvals,constants,options,numHeartBeats,indexes,simtime);
end
warning 'on'
disp('EOA simulated')

%% simulate range of ALVOT values
arange = linspace(minparamvalues(strcmp(paramNames,'A_ao')),maxparamvalues(strcmp(paramNames,'A_ao')) ,numsims);
arange = [minparamvalues(strcmp(paramNames,'A_ao')) - extrarange*minparamvalues(strcmp(paramNames,'A_ao')),...
    arange,...
    maxparamvalues(strcmp(paramNames,'A_ao')) + extrarange*maxparamvalues(strcmp(paramNames,'A_ao'))];

numHeartBeats = 20; %number of heartbeats to simulate
warning 'off'
pvals = parameterValuesOriginal;
asims = cell(size(arange));
for i = 1:length(arange)
    pvals(strcmp(paramNames,'A_ao')) = arange(i);
    [~,asims{i}] = simulate_avatarbelen(pvals,constants,options,numHeartBeats,indexes,simtime);
end
warning 'on'

disp('ELCo simulated')

%% simulate range of Emax_LV values
emaxrange = linspace(minparamvalues(strcmp(paramNames,'Emax_LV')),maxparamvalues(strcmp(paramNames,'Emax_LV')) ,numsims);
emaxrange = [minparamvalues(strcmp(paramNames,'Emax_LV')) - extrarange*minparamvalues(strcmp(paramNames,'Emax_LV')),...
    emaxrange,...
    maxparamvalues(strcmp(paramNames,'Emax_LV')) + extrarange*maxparamvalues(strcmp(paramNames,'Emax_LV'))];
warning 'off'
pvals = parameterValuesOriginal;
emaxsims = cell(size(emaxrange));
for i = 1:length(emaxrange)
    pvals(strcmp(paramNames,'Emax_LV')) = emaxrange(i);
    [~,emaxsims{i}] = simulate_avatarbelen(pvals,constants,options,numHeartBeats,indexes,simtime);
end
warning 'on'

disp('Emax LV simulated')

%% Create a results plot of the simulations with the range of parameter values
% define the names in the model
ynames = {'P_Aortic','Pperipheral','pressGrad MV','P Dmv','mv open','P Dav','av open','pressgrad AV','Ela','Elv','Qcaa','Qpc','P pulmvein','Qpvc','qLA','qLV','pLA','pLV','aaCorr','avCorr','mvCorr','pvCorr','Vla','Vlv','P_Brachial'};
xnames = {'Ppvc','Qpv','Vla','Qmv','Vlv','Qav','Paa','Qaa','Ppc','mv_open','av_open'};

% create the figure
figName='sensitivity_analysis';
fig=figure('Name',figName);
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 15;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])

sims = {eoasims,asims,emaxsims};
simnames = {'EOA','A_L_V_O_T','Emax LV'};
ranges = {eoarange,arange,emaxrange};
tiledlayout(length(sims),3,'TileSpacing','Compact','Padding','compact')


% create a default color map 
colorMapLength = length(eoarange);
darkblue = [0, 0, 255]./255;
skyeblue = [224, 240, 255]/255;
col1=skyeblue.*0.9;
col2=darkblue.*0.5;
colors_p = [linspace(col1(1),col2(1),colorMapLength)', linspace(col1(2),col2(2),colorMapLength)', linspace(col1(3),col2(3),colorMapLength)'];

% plot the figure
for s = 1:length(sims)
    sim = sims{s};
    nexttile
    hold on
    for i = 1:length(sim)
        pl(i) = plot(sim{i}.t,sim{i}.x(:,strcmp(xnames,'Qmv')),'-','LineWidth',1.1,'color',colors_p(i,:));
    end
    if s == 1
        title('Blood flow in the mitral valve')
    end
    xlabel('Time (s)')
    ylabel(sprintf('Sensitivity to %s\n\nBlood flow in\nmitral valve (ml/s)',simnames{s}))
    ylim([0,615])

    nexttile
    hold on
    for i = 1:length(sim)
        pl(i) = plot(sim{i}.t,sim{i}.x(:,strcmp(xnames,'Qav')),'-','LineWidth',1.1,'color',colors_p(i,:));
    end
    if s == 1
        title('Blood flow in the aortic valve')
    end
    xlabel('Time (s)')
    ylabel(sprintf('Blood flow in\naortic valve (ml/s)'))
    ylim([0,820])

    nexttile
    hold on
    for i = 1:length(sim)
        pl(i) = plot(sim{i}.t,sim{i}.y(:,strcmp(ynames,'P_Aortic')),'-','LineWidth',1.1,'color',colors_p(i,:));
    end
    if s == 1
        title('Aortic pressure')
    end
    xlabel('Time (s)')
    ylabel(sprintf('Aortic pressure\n(mmHg)'))
    ylim([40,150])
    legend([pl(1),pl(2),pl(end-1),pl(end)],{sprintf('%s = %0.2f (min-20%%)',simnames{s},ranges{s}(1)),sprintf('%s = %0.2f (min)',simnames{s},ranges{s}(2)),sprintf('%s = %0.2f (max)',simnames{s},ranges{s}(end-1)),sprintf('%s = %0.2f (max+20%%)',simnames{s},ranges{s}(end))},...
        'FontSize',7)
end

% Save the figure
exportgraphics(fig,[fullfile('Results',figName) '.pdf'],'ContentType','vector')
