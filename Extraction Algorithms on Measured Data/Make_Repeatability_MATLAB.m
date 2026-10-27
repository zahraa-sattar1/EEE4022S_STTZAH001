%% Repeated-walk error variation
% Zahraa Sattar | EEE4022S
%
% Plot the final results with sessions pooled by participant and pace.
% The P2_T9 exclusion is a separate diagnostic comparison.

root=fileparts(mfilename('fullpath'));
out=fullfile(root,'Report Results','Repeatability');
if ~isfolder(out), mkdir(out); end
T=readtable(fullfile(root,'Report Results','Results_AllConfigs.csv'),'TextType','string');
T=T(T.config=="Final" & T.walkType=="Single pass",:);
refs=[T.walkedDistance_m./T.referenceDuration_WAV_s,2*T.referenceDuration_WAV_s./T.gtStepCount,2*T.walkedDistance_m./T.gtStepCount];
E=100*([T.estSpeed_mps,T.estStrideDuration_mean_s,T.estStrideLength_mean_m]-refs)./refs;
people=["P1","P2","P3"]; paces=["Slow","Comfortable","Fast","Very fast"];
C=[.196 .424 .620;.725 .404 .141;.145 .502 .435];
S=zeros(3,4,3); SS=zeros(3,3); dfs=zeros(3,1); Swithout=zeros(4,3); ssWithout=zeros(1,3); dfWithout=0;
f=figure('Color','w','Units','centimeters','Position',[2 2 22 8]);
tl=tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
titles={'(a) Walking speed','(b) Stride duration','(c) Stride length'};
for j=1:3
 ax=nexttile; hold on;
 for p=1:3
  for k=1:4
   mask=T.participant==people(p) & T.pace==paces(k); vals=E(mask,j); centred=vals-mean(vals);
   S(p,k,j)=std(vals); SS(p,j)=SS(p,j)+sum(centred.^2);
   if j==1, dfs(p)=dfs(p)+numel(vals)-1; end
   jitter=linspace(-.065,.065,numel(vals))';
   scatter(k+(p-2)*.16+jitter,abs(centred),5,C(p,:),'filled','MarkerFaceAlpha',.28,'HandleVisibility','off');
  end
  plot((1:4)+(p-2)*.16,squeeze(S(p,:,j)),'-o','Color',C(p,:),'LineWidth',1.15,'MarkerSize',3,'DisplayName',people(p));
 end
 for k=1:4
  mask=T.participant=="P2" & T.pace==paces(k) & T.trial_label~="T9"; vals=E(mask,j);
  Swithout(k,j)=std(vals); ssWithout(j)=ssWithout(j)+sum((vals-mean(vals)).^2);
  if j==1,dfWithout=dfWithout+numel(vals)-1;end
 end
 plot(1:4,Swithout(:,j),'--s','Color',[.35 .18 .05],'LineWidth',1.15,'MarkerSize',3,'MarkerFaceColor','w','DisplayName','P2 without T9');
 title(titles{j},'FontSize',9,'FontWeight','bold');
 xlim([.65 4.35]);xticks(1:4);xticklabels({'Slow','Comf.','Fast','V. fast'});xtickangle(35);
 ylim([0 40]);yticks(0:10:40);grid on;ax.GridAlpha=.12;
 if j==1,ylabel({'Error variation','(percentage points)'});end
 xlabel('Pace');
 if j==3
  mask=T.participant=="P2" & T.pace=="Very fast"; v=E(mask,j); trial=T.trial_label(mask); yy=abs(v(trial=="T9")-mean(v));
  text(3.8,yy-4,'P2\_T9','HorizontalAlignment','right','FontSize',8,'Color',C(2,:));
 end
end
ax=nexttile;hold on;
pooled=sqrt(SS./dfs);
for p=1:3
 plot(1:3,pooled(p,:),'-o','Color',C(p,:),'LineWidth',1.2,'MarkerSize',4);
end
without=sqrt(ssWithout/dfWithout);
plot(1:3,without,'--s','Color',[.35 .18 .05],'LineWidth',1.15,'MarkerSize',4,'MarkerFaceColor','w');
title({'(d) Participant comparison','Within-pace pooled SD'},'FontSize',9,'FontWeight','bold');
xticks(1:3);xticklabels({'Speed','Duration','Length'});xtickangle(35);xlim([.7 3.3]);ylim([0 8]);yticks(0:2:8);ylabel('SD (percentage points)');grid on;ax.GridAlpha=.12;
ax1=findall(f,'Type','axes');
for a=ax1'
 set(a,'FontName','Times New Roman','FontSize',9,'Color','w','XColor','k','YColor','k','Box','off');axtoolbar(a,{});
end
axes(ax1(end));
h=gobjects(4,1);
for p=1:3,h(p)=plot(nan,nan,'-o','Color',C(p,:),'LineWidth',1.2,'MarkerSize',3);end
h(4)=plot(nan,nan,'--s','Color',[.35 .18 .05],'LineWidth',1.2,'MarkerSize',3);
lg=legend(h,{'P1','P2','P3','P2 without P2\_T9'},'Orientation','horizontal','Box','off','FontSize',9,'TextColor','k');lg.Layout.Tile='north';
xlabel(tl,{'Small dots: each walk''s absolute deviation from its group mean error.','Lines: SD. Dashed line: diagnostic exclusion of P2\_T9.'},'FontName','Times New Roman','FontSize',9,'Color','k');
set(findall(f,'Type','text'),'FontName','Times New Roman','Color','k');
exportgraphics(f,fullfile(out,'Fig_Measured_Repeatability_MATLAB.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(f,fullfile(out,'Fig_Measured_Repeatability_MATLAB.png'),'Resolution',300,'BackgroundColor','white');
savefig(f,fullfile(out,'Fig_Measured_Repeatability_MATLAB.fig'));
disp(pooled);disp(without);disp(Swithout);

pooledAll=sqrt(sum(SS,1)/sum(dfs));
writematrix(pooledAll,fullfile(out,'Pooled_SD.csv'));
writematrix(pooled,fullfile(out,'Participant_Pooled_SD.csv'));
