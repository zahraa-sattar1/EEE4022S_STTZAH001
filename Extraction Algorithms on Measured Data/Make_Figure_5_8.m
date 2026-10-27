%% Matched step-detection examples
% Zahraa Sattar | EEE4022S
%
% Redraw Figure 5.8 from the saved pipeline outputs for P1_T47 and P2_T9.
% This uses stored spectrograms and detections without reprocessing the recordings.

function Make_Figure_5_8
root=fileparts(mfilename('fullpath'));
source=fullfile(root,'Report Results','Undercount_Comparison.mat');
out=fullfile(root,'Report Results','Undercount Figure');
if ~isfolder(out), mkdir(out); end
load(source,'results','ids');
try, settings().matlab.appearance.figure.GraphicsTheme.TemporaryValue='light'; catch, end
f=figure('Color','w','Units','centimeters','Position',[2 2 32 13]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
summary=struct([]);
for k=1:2
 ax=nexttile(tl); p=results{k}.plot;
 t=p.legTime(:);v=p.F_smooth(:)*p.lambda/2;d=results{k}.stepTimes_s(:);
 imagesc(ax,t,p.legFreq*p.lambda/2,p.legPdB);axis(ax,'xy');hold(ax,'on');clim(ax,[-50 0]);colormap(ax,turbo);
 plot(ax,t,v,'w','LineWidth',1.6);
 scatter(ax,d,interp1(t,v,d),55,'o','MarkerEdgeColor',[1 .25 .25],'LineWidth',1.5);
 ylim(ax,[0 5]);xlim(ax,[0 3.6]);
 set(ax,'FontSize',16,'FontName','Times New Roman','Color','w','XColor','k','YColor','k');
 xlabel(ax,'Time (s)','FontSize',16);ylabel(ax,'Radial speed magnitude (m/s)','FontSize',16);
 title(ax,sprintf('(%c) %s: reference 1.41m/s, %d / 8','a'+k-1,strrep(ids{k},'_','\_'),results{k}.stepCount),'Color','k','FontSize',15);
 axtoolbar(ax,{});
 callout(ax,.45,4.70,.67,3.45,'Leg return');
 callout(ax,.50,.40,.85,1.30,'Torso return');
 callout(ax,1.68,.40,1.30,.70,{'Arm return','(likely)'});
 peakV=interp1(t,v,d);depth=zeros(numel(d)-1,1);
 for j=1:numel(d)-1
  depth(j)=min(peakV(j:j+1))-min(v(t>=d(j)&t<=d(j+1)));
 end
 summary(k).id=ids{k};summary(k).stepTimes_s=d';summary(k).dipDepths_mps=depth';
 summary(k).medianDipDepth_mps=median(depth);summary(k).largestGap_s=max(diff(d));
 if k==1
  label(ax,2.78,4.65,sprintf('Median dip depth\n%.2f m/s',median(depth)));
 else
  label(ax,2.1,4.32,'Poorly separated swings');
  label(ax,2.91,3.73,sprintf('Median dip depth\n%.1f m/s',median(depth)));
  x1=d(3);x2=d(4);y=4.65;
  plot(ax,[x1 x1],[interp1(t,v,x1)+.08 y],'w:','LineWidth',1);
  plot(ax,[x2 x2],[interp1(t,v,x2)+.08 y],'w:','LineWidth',1);
  doubleArrow(ax,x1,x2,y);
  label(ax,(x1+x2)/2,4.86,sprintf('Gap: %.3f s',x2-x1));
  summary(k).annotatedCentralGap_s=x2-x1;
 end
end
cb=colorbar(ax);cb.Label.String={'Power relative to guarded','reference (dB)'};
cb.FontSize=15;cb.Label.FontSize=15;cb.Color='k';cb.Label.Color='k';
h1=plot(ax,nan,nan,'k','LineWidth',1.6);h2=plot(ax,nan,nan,'o','Color',[1 .25 .25],'LineWidth',1.5);
lg=legend(ax,[h1 h2],'Processed leg-speed curve','Final detections','Orientation','horizontal','FontSize',16);
lg.Layout.Tile='south';lg.TextColor='k';lg.Color='w';
set(findall(f,'Type','text'),'FontName','Times New Roman');drawnow;
name='Fig_StepCount_MatchedComparison_v4';
savefig(f,fullfile(out,[name '.fig']));
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',300,'BackgroundColor','white');
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
fid=fopen(fullfile(out,'Annotation_Checks.json'),'w');fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));fclose(fid);
close(f);
end
function label(ax,x,y,str)
text(ax,x,y,str,'Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[.09 .06 .15],'Margin',2);
end
function callout(ax,x,y,tx,ty,str)
dx=tx-x;dy=ty-y;quiver(ax,x,y,dx,dy,0,'Color','w','LineWidth',1.3,'MaxHeadSize',.25,'HandleVisibility','off');
label(ax,x,y,str);
end
function doubleArrow(ax,x1,x2,y)
plot(ax,[x1 x2],[y y],'w','LineWidth',1.3);
plot(ax,[x1+.055 x1 x1+.055],[y+.055 y y-.055],'w','LineWidth',1.3);
plot(ax,[x2-.055 x2 x2-.055],[y+.055 y y-.055],'w','LineWidth',1.3);
end
