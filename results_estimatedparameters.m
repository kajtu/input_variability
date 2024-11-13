clear
close all

addpath('Requirements')
addpath('Data')

%% Load data 
% (NOTE: data not provided in gitgub due to ethical restrictions)
p_intrainter = readtable('Analysis_results_BCvsFV.xlsx','Sheet','Parameter results_intra_inter');
paramnames = p_intrainter.Var1(1:23);

p_intrainter_m = readmatrix('Analysis_results_BCvsFV.xlsx','Sheet','Parameter results_intra_inter');
p_intrainter_m = p_intrainter_m(:,4:33);
p_bc1 = p_intrainter_m(:,1:3:30);
p_bc2 = p_intrainter_m(:,2:3:30);
p_fv = p_intrainter_m(:,3:3:30);
p_bc1=p_bc1';
p_bc2=p_bc2';
p_fv=p_fv';

p_seq_m= readmatrix('Analysis_results_BCvsFV.xlsx','Sheet','Parameter results EPI TFE');
p_seq_m = p_seq_m(:,4:23);
p_epi =p_seq_m(:,2:2:end);
p_tfe = p_seq_m(:,1:2:end);
p_epi = p_epi';
p_tfe=p_tfe';

%% allocate variables
cov_intra=nan(size(paramnames));
cov_inter=nan(size(paramnames));
cov_seq=nan(size(paramnames));
meanDiff_intra=nan(size(paramnames));
meanp196D_intra=nan(size(paramnames));
meanm196D_intra=nan(size(paramnames));
mean_intra=nan(size(paramnames));
sd_intra=nan(size(paramnames));
p_intra_t=nan(size(paramnames));
p_intra_mw=nan(size(paramnames));
meanDiffPercent_intra =nan(size(paramnames));

meanDiff_inter=nan(size(paramnames));
meanp196D_inter=nan(size(paramnames));
meanm196D_inter=nan(size(paramnames));
mean_inter=nan(size(paramnames));
sd_inter=nan(size(paramnames));
p_inter_t=nan(size(paramnames));
p_inter_mw=nan(size(paramnames));
meanDiffPercent_inter =nan(size(paramnames));


meanDiff_seq=nan(size(paramnames));
meanp196D_seq=nan(size(paramnames));
meanm196D_seq=nan(size(paramnames));
mean_seq=nan(size(paramnames));
sd_seq=nan(size(paramnames));
p_seq_t=nan(size(paramnames));
p_seq_mw=nan(size(paramnames));
meanDiffPercent_seq =nan(size(paramnames));

m_intra=cell(size(paramnames));
interval_intra=cell(size(paramnames));
intra95percent=cell(size(paramnames));
m_inter=cell(size(paramnames));
interval_inter=cell(size(paramnames));
inter95percent=cell(size(paramnames));
m_seq=cell(size(paramnames));
interval_seq=cell(size(paramnames));
seq95percent=cell(size(paramnames));


cov_intra_t=cell(size(paramnames));
cov_inter_t=cell(size(paramnames));
cov_seq_t=cell(size(paramnames));

p_intra=cell(size(paramnames));
p_inter=cell(size(paramnames));
p_seq=cell(size(paramnames));


%% For each parameter, compute bland-altman analysis, COV, and p-values
string1='mean (analysis1 & analysis2)';
string2='analysis1 - analysis2';
for p = 1:length(paramnames)
    cov_intra(p) = compute_COV(p_bc1(:,p),p_bc2(:,p));
    cov_inter(p) = compute_COV(p_bc1(:,p),p_fv(:,p));
    cov_seq(p) = compute_COV(p_epi(:,p),p_tfe(:,p));
    cov_intra_t{p} = sprintf('%.2g',cov_intra(p));
    cov_inter_t{p} = sprintf('%.2g',cov_inter(p));
    cov_seq_t{p} = sprintf('%.2g',cov_seq(p));

  
    [meanDiff_intra(p),meanp196D_intra(p),meanm196D_intra(p),meanDiffPercent_intra(p)]=BlandAltmanpar(p_bc1(:,p),p_bc2(:,p),string1,string2);
    title(paramnames{p})
    mean_intra(p) = mean([p_bc1(:,p);p_bc2(:,p)]);
    sd_intra(p) = std([p_bc1(:,p);p_bc2(:,p)]);
    [~,p_intra_t(p)] = ttest2(p_bc1(:,p),p_bc2(:,p));
    p_intra_mw(p) = ranksum(p_bc1(:,p),p_bc2(:,p));%,'method','exact'
    m_intra{p} = sprintf('%.2g +- %.2g',mean_intra(p),sd_intra(p));
    intra95percent{p} = sprintf('%.2g\n(%.2g to %.2g)',meanDiff_intra(p),meanm196D_intra(p)/sqrt(10),meanp196D_intra(p)/sqrt(10));
    interval_intra{p} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_intra(p),meanm196D_intra(p),meanp196D_intra(p));
    p_intra{p} = sprintf('%.2g',p_intra_mw(p));

    [meanDiff_inter(p),meanp196D_inter(p),meanm196D_inter(p),meanDiffPercent_inter(p)]=BlandAltmanpar(p_bc1(:,p),p_fv(:,p),string1,string2);
    title(paramnames{p})
    mean_inter(p) = mean([p_bc1(:,p);p_fv(:,p)]);
    sd_inter(p) = std([p_bc1(:,p);p_fv(:,p)]);
    [~,p_inter_t(p)] = ttest2(p_bc1(:,p),p_fv(:,p));
    p_inter_mw(p) = ranksum(p_bc1(:,p),p_fv(:,p));
    m_inter{p} = sprintf('%.2g +- %.2g',mean_inter(p),sd_inter(p));
    inter95percent{p} = sprintf('%.2g\n(%.2g to %.2g)',meanDiff_inter(p),meanm196D_inter(p)/sqrt(10),meanp196D_inter(p)/sqrt(10));
    interval_inter{p} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_inter(p),meanm196D_inter(p),meanp196D_inter(p));
    p_inter{p} = sprintf('%.2g',p_inter_mw(p));

    [meanDiff_seq(p),meanp196D_seq(p),meanm196D_seq(p),meanDiffPercent_seq(p)]=BlandAltmanpar(p_tfe(:,p),p_epi(:,p),string1,string2);
    title(paramnames{p})
    mean_seq(p) = mean([p_epi(:,p);p_tfe(:,p)]);
    sd_seq(p) = std([p_epi(:,p);p_tfe(:,p)]);
    [~,p_seq_t(p)] = ttest2(p_epi(:,p),p_tfe(:,p));
    p_seq_mw(p) = ranksum(p_epi(:,p),p_tfe(:,p));
    m_seq{p} = sprintf('%.2g +- %.2g',mean_seq(p),sd_seq(p));
    seq95percent{p} = sprintf('%.2g\n(%.2g to %.2g)',meanDiff_seq(p),meanm196D_seq(p)/sqrt(10),meanp196D_seq(p)/sqrt(10));
    interval_seq{p} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_seq(p),meanm196D_seq(p),meanp196D_seq(p));
    p_seq{p} = sprintf('%.2g',p_seq_mw(p));

    close all
end

%% Do the same for all selected parameters
selectedpars = {'Emin_LV','Emax_LV','Caa','m1_LV','m2_LV','k_diast_LV','k_syst_LV'};
for p = 1:length(selectedpars)
    indp = strcmp(selectedpars{p},paramnames);
    string1='mean (analysis1 & analysis2)';
    string2='analysis1 - analysis2';
    BlandAltmanpar(p_bc1(:,indp),p_bc2(:,indp),string1,string2);
    title(selectedpars{p})
    string1='mean (observer1 & observer2)';
    string2='observer1 - observer2';
    BlandAltmanpar(p_bc1(:,indp),p_fv(:,indp),string1,string2);
    title(selectedpars{p})
    string1='mean (SGRE & EPI)';
    string2='SGRE - EPI';
    BlandAltmanpar(p_tfe(:,indp),p_epi(:,indp),string1,string2);
    title(selectedpars{p})
    if strcmp(selectedpars{p},'Emax_LV')
        ylim([-0.8,0.8])
        yticks([-0.8:0.2:0.8])
        axis('square')
        fig = figure(6);
        fontsize(fig,12,'points')
    end
end
close all

%% Print out relevant statistics
inds=ismember(paramnames,{'Emin_LV','Emax_LV','Caa','m1_LV','m2_LV','k_diast_LV','k_syst_LV'});
inds = find(inds);

fprintf('CoV for interobserver is larger than intraobserver variability for all the 7 selected params?: %d\n',sum(cov_inter(inds)>cov_intra(inds)) == length(inds))
disp('-------------')

fprintf('Smallest p-value for ALL params intraobserver variability: %d\n',min(p_intra_mw))
fprintf('Smallest p-value for 7 selected params intraobserver variability: %d\n',min(p_intra_mw(inds)))
fprintf('Smallest p-value for ALL params interobserver variability: %d\n',min(p_inter_mw))
fprintf('Smallest p-value for 7 selected params interobserver variability: %d\n',min(p_inter_mw(inds)))
fprintf('Smallest p-value for ALL params intersequence variability: %d\n',min(p_seq_mw))
fprintf('Smallest p-value for 7 selected params intersequence variability: %d\n',min(p_seq_mw(inds)))
disp('-------------')

selectedparamnames =paramnames(inds);
[m,ind]=max(cov_inter);
fprintf('Largest CoV for interobserver among ALL params: %s (%0.2f %%)\n',paramnames{ind},m)
[m,ind]=max(cov_inter(inds));
fprintf('Largest CoV for interobserver among the 7 selected params: %s (%0.2f %%)\n',selectedparamnames{ind},m)
[m,ind]=min(cov_inter);
fprintf('Smallest CoV for interobserver among ALL params: %s (%0.2f %%)\n',paramnames{ind},m)
[m,ind]=min(cov_inter(inds));
fprintf('Smallest CoV for interobserver among the 7 selected params: %s (%0.2f %%)\n',selectedparamnames{ind},m)
disp('-------------')
[m,ind]=max(cov_intra);
fprintf('Largest CoV for intraobserver among ALL params: %s (%0.2f %%)\n',paramnames{ind},m)
[m,ind]=max(cov_intra(inds));
fprintf('Largest CoV for intraobserver among the 7 selected params: %s (%0.2f %%)\n',selectedparamnames{ind},m)
[m,ind]=min(cov_intra);
fprintf('Smallest CoV for intraobserver among ALL params: %s (%0.2f %%)\n',paramnames{ind},m)
[m,ind]=min(cov_intra(inds));
fprintf('Smallest CoV for intraobserver among the 7 selected params: %s (%0.2f %%)\n',selectedparamnames{ind},m)
disp('-------------')
[m,ind]=max(cov_seq);
fprintf('Largest CoV for intersequence among ALL params: %s (%0.2f %%)\n',paramnames{ind},m)
[m,ind]=max(cov_seq(inds));
fprintf('Largest CoV for intersequence among the 7 selected params: %s (%0.2f %%)\n',selectedparamnames{ind},m)
[m,ind]=min(cov_seq);
fprintf('Smallest CoV for intersequence among ALL params: %s (%0.2f %%)\n',paramnames{ind},m)
[m,ind]=min(cov_seq(inds));
fprintf('Smallest CoV for intersequence among the 7 selected params: %s (%0.2f %%)\n',selectedparamnames{ind},m)

%% Write tables to excel files
folder = 'Results';
varnames = {'mean +- SD','mean difference (%)','mean difference (-1.96*SD,+1.96*SD)','p-value','CoV (%)'};
intra=table(m_intra,meanDiffPercent_intra,interval_intra,p_intra,cov_intra_t,'Rownames',paramnames,'Variablenames',varnames);
inter=table(m_inter,meanDiffPercent_inter,interval_inter,p_inter,cov_inter_t,'Rownames',paramnames,'Variablenames',varnames);
seq=table(m_seq,interval_seq,meanDiffPercent_seq,p_seq,cov_seq_t,'Rownames',paramnames,'Variablenames',varnames);

writetable(intra,fullfile(folder,'params_all_intra.xlsx'),"WriteRowNames",1)
writetable(inter,fullfile(folder,'params_all_inter.xlsx'),"WriteRowNames",1)
writetable(seq,fullfile(folder,'params_all_seq.xlsx'),"WriteRowNames",1)

writetable(intra(inds,:),fullfile(folder,'params_selected_intra.xlsx'),"WriteRowNames",1)
writetable(inter(inds,:),fullfile(folder,'params_selected_inter.xlsx'),"WriteRowNames",1)
writetable(seq(inds,:),fullfile(folder,'params_selected_seq.xlsx'),"WriteRowNames",1)

disp('Done. Resulting tables are stored in the Results folder.')
