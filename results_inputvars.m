clear
close all

addpath('Data')
addpath('Requirements')

%% Load data
% (NOTE: data not provided in gitgub due to ethical restrictions)
input = readmatrix('Analysis_results_BCvsFV.xlsx','Sheet','Comparison_Aao_VTI_ESV_SV');
folder='Results';
 
input = input(1:12,2:33);
paramnames = {'ALVOT','ESV','T','EOA'}';
input = input(~isnan(input(:,1)),:);
p_bc1 = input(:,[1,9,14,17]);
p_bc2 = input(:,[2,10,14,18]);
p_fv = input(:,[3,11,14,19]);
p_tfe = input(:,[3,11,12,19]);
p_epi = input(:,[4,14,13,20]);

%% Allocate variables 
cov_intra=nan(size(paramnames));
cov_inter=nan(size(paramnames));
cov_seq=nan(size(paramnames));
meanDiff_intra=nan(size(paramnames));
meanp196D_intra=nan(size(paramnames));
meanm196D_intra=nan(size(paramnames));
mean_intra=nan(size(paramnames));
sd_intra=nan(size(paramnames));
p_intra_mw=nan(size(paramnames));

meanDiff_inter=nan(size(paramnames));
meanp196D_inter=nan(size(paramnames));
meanm196D_inter=nan(size(paramnames));
mean_inter=nan(size(paramnames));
sd_inter=nan(size(paramnames));
p_inter_mw=nan(size(paramnames));

meanDiff_seq=nan(size(paramnames));
meanp196D_seq=nan(size(paramnames));
meanm196D_seq=nan(size(paramnames));
mean_seq=nan(size(paramnames));
sd_seq=nan(size(paramnames));
p_seq_mw=nan(size(paramnames));

m_intra=cell(size(paramnames));
interval_intra=cell(size(paramnames));
intra95percent=cell(size(paramnames));
m_inter=cell(size(paramnames));
interval_inter=cell(size(paramnames));
inter95percent=cell(size(paramnames));
m_seq=cell(size(paramnames));
interval_seq=cell(size(paramnames));

cov_intra_t=cell(size(paramnames));
cov_inter_t=cell(size(paramnames));
cov_seq_t=cell(size(paramnames));

p_intra=cell(size(paramnames));
p_inter=cell(size(paramnames));
p_seq=cell(size(paramnames));

%% For each input variable, compute bland-altman analysis, COV, and p-values
string1='mean (analysis1 & analysis2)';
string2='analysis1 - analysis2';
for p = 1:length(paramnames)
    cov_intra(p) = compute_COV(p_bc1(:,p),p_bc2(:,p));
    cov_inter(p) = compute_COV(p_bc1(:,p),p_fv(:,p));
    cov_seq(p) = compute_COV(p_epi(:,p),p_tfe(:,p));
    cov_intra_t{p} = sprintf('%.2g',cov_intra(p));
    cov_inter_t{p} = sprintf('%.2g',cov_inter(p));
    cov_seq_t{p} = sprintf('%.2g',cov_seq(p));

    if sum(isnan(p_bc1(:,p)))<1
        [meanDiff_intra(p),meanp196D_intra(p),meanm196D_intra(p)]=BlandAltmanpar(p_bc1(:,p),p_bc2(:,p),string1,string2);
        title(paramnames{p})
        mean_intra(p) = mean([p_bc1(:,p);p_bc2(:,p)]);
        sd_intra(p) = std([p_bc1(:,p);p_bc2(:,p)]);
        m_intra{p} = sprintf('%.2g +- %.2g',mean_intra(p),sd_intra(p));
        interval_intra{p} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_intra(p),meanm196D_intra(p),meanp196D_intra(p));
        p_intra{p} = sprintf('%.2g',p_intra_mw(p));

        [meanDiff_inter(p),meanp196D_inter(p),meanm196D_inter(p)]=BlandAltmanpar(p_bc1(:,p),p_fv(:,p),string1,string2);
        title(paramnames{p})
        mean_inter(p) = mean([p_bc1(:,p);p_fv(:,p)]);
        sd_inter(p) = std([p_bc1(:,p);p_fv(:,p)]);
        p_inter_mw(p) = ranksum(p_bc1(:,p),p_fv(:,p));
        m_inter{p} = sprintf('%.2g +- %.2g',mean_inter(p),sd_inter(p));
        interval_inter{p} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_inter(p),meanm196D_inter(p),meanp196D_inter(p));
        p_inter{p} = sprintf('%.2g',p_inter_mw(p));
    end

    if sum(isnan(p_epi(:,p)))<1
        [meanDiff_seq(p),meanp196D_seq(p),meanm196D_seq(p)]=BlandAltmanpar(p_tfe(:,p),p_epi(:,p),string1,string2);
        title(paramnames{p})
        mean_seq(p) = mean([p_epi(:,p);p_tfe(:,p)]);
        sd_seq(p) = std([p_epi(:,p);p_tfe(:,p)]);
        p_seq_mw(p) = ranksum(p_epi(:,p),p_tfe(:,p));
        m_seq{p} = sprintf('%.2g +- %.2g',mean_seq(p),sd_seq(p));
        interval_seq{p} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_seq(p),meanm196D_seq(p),meanp196D_seq(p));
        p_seq{p} = sprintf('%.2g',p_seq_mw(p));
    end
    close all
end


%% Save tables in excel files
varnames = {'mean +- SD','mean difference (-1.96*SD,+1.96*SD)','p-value','CoV (%)'};
intra=table(m_intra,interval_intra,p_intra,cov_intra_t,'Rownames',paramnames,'Variablenames',varnames);
inter=table(m_inter,interval_inter,p_inter,cov_inter_t,'Rownames',paramnames,'Variablenames',varnames);
seq=table(m_seq,interval_seq,p_seq,cov_seq_t,'Rownames',paramnames,'Variablenames',varnames);

writetable(intra,fullfile(folder,'input_all_intra.xlsx'),"WriteRowNames",1)
writetable(inter,fullfile(folder,'input_all_inter.xlsx'),"WriteRowNames",1)
writetable(seq,fullfile(folder,'input_all_seq.xlsx'),"WriteRowNames",1)

%% Do analysis of total flows
% Load data
load('mv_totflow.mat')
load('av_totflow.mat')
load('aa_totflow.mat')

% Compute statistics for MV
cov_mv_intra = compute_COV(mv_totflow.BC11,mv_totflow.BC12);
cov_mv_inter = compute_COV(mv_totflow.BC11,mv_totflow.FV5);
cov_mv_seq = compute_COV(mv_totflow.FV5,mv_totflow.EPI1);
p_mv_mw(1) = ranksum(mv_totflow.BC11,mv_totflow.BC12);
p_mv_mw(2) = ranksum(mv_totflow.BC11,mv_totflow.FV5);
p_mv_mw(3) = ranksum(mv_totflow.FV5,mv_totflow.EPI1);
[meanDiff_mv(1),meanp196D_mv(1),meanm196D_mv(1),meanDiffPercent_mv(1)]=BlandAltmanpar(mv_totflow.BC11,mv_totflow.BC12,'mean of analysis 1 & 2','Analysis 1 - Analysis 2');
title('MV net flow volume')
[meanDiff_mv(2),meanp196D_mv(2),meanm196D_mv(2),meanDiffPercent_mv(2)]=BlandAltmanpar(mv_totflow.BC11,mv_totflow.FV5,'mean of observer 1 & 2','Observer 1 - Observer 2');
title('MV net flow volume')
[meanDiff_mv(3),meanp196D_mv(3),meanm196D_mv(3),meanDiffPercent_mv(3)]=BlandAltmanpar(mv_totflow.FV5,mv_totflow.EPI1,'mean of SGRE & EPI','SGRE-EPI');
title('MV net flow volume')
m_intra_flow{1} = sprintf('%.2g +- %.2g',mean([mv_totflow.BC11;mv_totflow.BC12]),std([mv_totflow.BC11;mv_totflow.BC12]));
m_inter_flow{1} = sprintf('%.2g +- %.2g',mean([mv_totflow.BC11;mv_totflow.FV5]),std([mv_totflow.BC11;mv_totflow.FV5]));
m_seq_flow{1} = sprintf('%.2g +- %.2g',mean([mv_totflow.FV5;mv_totflow.EPI1]),std([mv_totflow.FV5;mv_totflow.EPI1]));

% Compute statistics for AV
cov_av_intra = compute_COV(av_totflow.BC13,av_totflow.BC14);
cov_av_inter = compute_COV(av_totflow.BC13,av_totflow.FV6);
cov_av_seq = compute_COV(av_totflow.FV6,av_totflow.EPI2);
p_av_mw(1) = ranksum(av_totflow.BC13,av_totflow.BC14);
p_av_mw(2) = ranksum(av_totflow.BC13,av_totflow.FV6);
p_av_mw(3) = ranksum(av_totflow.FV6,av_totflow.EPI2);
[meanDiff_av(1),meanp196D_av(1),meanm196D_av(1),meanDiffPercent_av(1)]=BlandAltmanpar(av_totflow.BC13,av_totflow.BC14,'mean of analysis 1 & 2','Analysis 1 - Analysis 2');
title('AV net flow volume')
[meanDiff_av(2),meanp196D_av(2),meanm196D_av(2),meanDiffPercent_av(2)]=BlandAltmanpar(av_totflow.BC13,av_totflow.FV6,'mean of observer 1 & 2','Observer 1 - Observer 2');
title('AV net flow volume')
[meanDiff_av(3),meanp196D_av(3),meanm196D_av(3),meanDiffPercent_av(3)]=BlandAltmanpar(av_totflow.FV6,av_totflow.EPI2,'mean of SGRE & EPI','SGRE-EPI');
title('AV net flow volume')
m_intra_flow{2} = sprintf('%.2g +- %.2g',mean([av_totflow.BC13;av_totflow.BC14]),std([av_totflow.BC13;av_totflow.BC14]));
m_inter_flow{2} = sprintf('%.2g +- %.2g',mean([av_totflow.BC13;av_totflow.FV6]),std([av_totflow.BC13;av_totflow.FV6]));
m_seq_flow{2} = sprintf('%.2g +- %.2g',mean([av_totflow.FV6;av_totflow.EPI2]),std([av_totflow.FV6;av_totflow.EPI2]));

% Compute statistics for AA
cov_aa_intra = compute_COV(aa_totflow.BC15,aa_totflow.BC16);
cov_aa_inter = compute_COV(aa_totflow.BC15,aa_totflow.FV7);
cov_aa_seq = compute_COV(aa_totflow.FV7,aa_totflow.EPI3);
p_aa_mw(1) = ranksum(aa_totflow.BC15,aa_totflow.BC16);
p_aa_mw(2) = ranksum(aa_totflow.BC15,aa_totflow.FV7);
p_aa_mw(3) = ranksum(aa_totflow.FV7,aa_totflow.EPI3);
[meanDiff_av(1),meanp196D_aa(1),meanm196D_aa(1),meanDiffPercent_aa(1)]=BlandAltmanpar(aa_totflow.BC15,aa_totflow.BC16,'mean of analysis 1 & 2','Analysis 1 - Analysis 2');
title('AA net flow volume')
[meanDiff_aa(2),meanp196D_aa(2),meanm196D_aa(2),meanDiffPercent_aa(2)]=BlandAltmanpar(aa_totflow.BC15,aa_totflow.FV7,'mean of observer 1 & 2','Observer 1 - Observer 2');
title('AA net flow volume')
[meanDiff_aa(3),meanp196D_aa(3),meanm196D_aa(3),meanDiffPercent_aa(3)]=BlandAltmanpar(aa_totflow.FV7,aa_totflow.EPI3,'mean of SGRE & EPI','SGRE-EPI');
title('AA net flow volume')
m_intra_flow{3} = sprintf('%.2g +- %.2g',mean([aa_totflow.BC15;aa_totflow.BC16]),std([aa_totflow.BC15;aa_totflow.BC16]));
m_inter_flow{3} = sprintf('%.2g +- %.2g',mean([aa_totflow.BC15;aa_totflow.FV7]),std([aa_totflow.BC15;aa_totflow.FV7]));
m_seq_flow{3} = sprintf('%.2g +- %.2g',mean([aa_totflow.FV7;aa_totflow.EPI3]),std([aa_totflow.FV7;aa_totflow.EPI3]));

close all

% Combine into one table and save to excel files
cov_table = table([cov_mv_intra,cov_mv_inter,cov_mv_seq]', ...
    [cov_av_intra,cov_av_inter,cov_av_seq]', ...
    [cov_aa_intra,cov_aa_inter,cov_aa_seq]','VariableNames',{'mv','av','aa'},'rownames',{'intra','inter','seq'});

pval_table = table(p_mv_mw',p_av_mw',p_aa_mw','VariableNames',{'mv','av','aa'},'rownames',{'intra','inter','seq'});

writetable(cov_table,fullfile(folder,'flows_cov_table.xlsx'),"WriteRowNames",1)
writetable(pval_table,fullfile(folder,'flows_pval_table.xlsx'),"WriteRowNames",1)


%% Create supplementary table of total flow statistics for revision 2
cov_intra_flow{1} = sprintf('%.2g',cov_mv_intra);
cov_intra_flow{2} = sprintf('%.2g',cov_av_intra);
cov_intra_flow{3} = sprintf('%.2g',cov_aa_intra);
cov_inter_flow{1} = sprintf('%.2g',cov_mv_inter);
cov_inter_flow{2} = sprintf('%.2g',cov_av_inter);
cov_inter_flow{3} = sprintf('%.2g',cov_aa_inter);
cov_seq_flow{1} = sprintf('%.2g',cov_mv_seq);
cov_seq_flow{2} = sprintf('%.2g',cov_av_seq);
cov_seq_flow{3} = sprintf('%.2g',cov_aa_seq);

%mv
interval_intra_flow{1} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_mv(1),meanm196D_mv(1),meanp196D_mv(1));
interval_inter_flow{1} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_mv(2),meanm196D_mv(2),meanp196D_mv(2));
interval_seq_flow{1} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_mv(3),meanm196D_mv(3),meanp196D_mv(3));
%av
interval_intra_flow{2} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_av(1),meanm196D_av(1),meanp196D_av(1));
interval_inter_flow{2} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_av(2),meanm196D_av(2),meanp196D_av(2));
interval_seq_flow{2} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_av(3),meanm196D_av(3),meanp196D_av(3));
%aa
interval_intra_flow{3} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_aa(1),meanm196D_aa(1),meanp196D_aa(1));
interval_inter_flow{3} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_aa(2),meanm196D_aa(2),meanp196D_aa(2));
interval_seq_flow{3} = sprintf('%.2g\n(%.2g, %.2g)',meanDiff_aa(3),meanm196D_aa(3),meanp196D_aa(3));

%mv
p_intra_flow{1} = sprintf('%.2g',p_mv_mw(1));
p_inter_flow{1} = sprintf('%.2g',p_mv_mw(2));
p_seq_flow{1} = sprintf('%.2g',p_mv_mw(3));
%av
p_intra_flow{2} = sprintf('%.2g',p_av_mw(1));
p_inter_flow{2} = sprintf('%.2g',p_av_mw(2));
p_seq_flow{2} = sprintf('%.2g',p_av_mw(3));
%aa
p_intra_flow{3} = sprintf('%.2g',p_aa_mw(1));
p_inter_flow{3} = sprintf('%.2g',p_aa_mw(2));
p_seq_flow{3} = sprintf('%.2g',p_aa_mw(3));

varnames = {'mean +- SD','mean difference (-1.96*SD,+1.96*SD)','p-value','CoV (%)'};
intra=table(m_intra_flow',interval_intra_flow',p_intra_flow',cov_intra_flow','Rownames',{'MV','AV','AA'},'Variablenames',varnames);
inter=table(m_inter_flow',interval_inter_flow',p_inter_flow',cov_inter_flow','Rownames',{'MV','AV','AA'},'Variablenames',varnames);
seq=table(m_seq_flow',interval_seq_flow',p_seq_flow',cov_seq_flow','Rownames',{'MV','AV','AA'},'Variablenames',varnames);

writetable(intra,fullfile(folder,'flows_all_intra.xlsx'),"WriteRowNames",1)
writetable(inter,fullfile(folder,'flows_all_inter.xlsx'),"WriteRowNames",1)
writetable(seq,fullfile(folder,'flows_all_seq.xlsx'),"WriteRowNames",1)

%% RMSE
%load data
input = readmatrix('Analysis_results_BCvsFV.xlsx','Sheet','RMSE');
mean = input(13,2:end);
sd = input(14,2:end);

% plot
meanvals = [mean(1:3:9);mean(2:3:9);mean(3:3:9)];
sdvals = [sd(1:3:9);sd(2:3:9);sd(3:3:9)];
X = categorical({'Mitral valve','Aortic valve','Ascending aorta'});
X = reordercats(X,{'Mitral valve','Aortic valve','Ascending aorta'});

fig=figure();
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 8;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
hold on
b=bar(X,meanvals,'linewidth',1);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0.5 0.5];
b(3).FaceColor = [1 1 1];
% colororder("tab")
a=[0.78,1,1.23];
errorbar([a;a+1;a+2],meanvals,sdvals,'k','linestyle','none','linewidth',1)
legend({'Intra-observer','Inter-observer','SGRE vs EPI'},'fontsize',12)
ylim([-5,85])
yticks(0:20:80)
ylabel('Root mean square error (mL/s)')
set(gca,'fontsize',12)

% Save the figure
exportgraphics(fig,[fullfile('Results','RMSE'), '.pdf'],'ContentType','vector')

%%
disp('Done. Resulting tables are stored in the Results folder.')
