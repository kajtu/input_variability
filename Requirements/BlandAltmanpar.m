function [meanDiff,meanp196D,meanm196D,meanDiffPercent] = BlandAltmanpar(A,B,string1,string2)
AB=[A;B];
meanAB=(A+B)./2; 
difff=A-B;
meanDiff=mean(difff);
stdDiff=std(difff);

meanDiffPercent = mean(100*difff./meanAB);

meanp196D=meanDiff+1.96*stdDiff;
meanm196D=meanDiff-1.96*stdDiff;
n=length(difff);
minD=min(meanAB)-0.1;
maxD=max(meanAB)+0.1;

figure

h1=line(meanAB,difff);
set(h1                         , ...
  'LineStyle'       , 'none'      , ...
  'Marker'          , '.'         , ...
  'MarkerSize'      , 7           , ...
  'Color'           , 'b'         , ...
  'LineWidth'       , 1        );

linelimit=max(AB);
linelimit2=min(AB);

scatter(meanAB,difff,'filled')
hold on;
plot([linelimit2; linelimit],ones(1,2)*meanp196D,'--r');
% hText1=text(linelimit-10,meanp196D+1.5,'+1.96 SD');
% hText2=text(linelimit-10,meanp196D+0.75, sprintf('%.3f',meanp196D));
hold on;
plot([linelimit2; linelimit],ones(1,2)*meanm196D,'--r');
% hText3=text(linelimit-10,meanm196D-0.75,'-1.96 SD');
% hText4=text(linelimit-10,meanm196D-1.5, sprintf('%.3f',meanm196D));
hold on;
plot([linelimit2; linelimit],ones(1,2)*meanDiff,'--','color',[0.5 0.5 0.5]);
% hText5=text(linelimit-10,meanDiff+0.75,'mean');
% hText6=text(linelimit-10,meanDiff-0.75, sprintf('%.3f',meanDiff));

plot([linelimit2; linelimit],zeros(size(ones(1,2)*meanDiff)),'k--')

xlim([linelimit2 linelimit]);
% ylim([-7 7]);
hXLabel=xlabel(string1);
hYLabel=ylabel(string2);

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 14   , ...
    'Fontweight' , 'bold');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'fontsize'    , 14, ...
  'LineWidth'   , 1         , ...
  'Fontweight' , 'bold');
