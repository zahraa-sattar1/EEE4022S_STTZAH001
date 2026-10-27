%% Measured development stages and speed benchmark
% Zahraa Sattar | EEE4022S
%
% Recalculate the report comparisons from the saved per-walk and per-trial results.
% The centroid benchmark uses the earlier run quoted in Table 5.5.

function Report_Comparisons
base=fullfile(fileparts(mfilename('fullpath')),'Report Results');
out=fullfile(base,'Comparisons');
if ~isfolder(out), mkdir(out); end
T=readtable(fullfile(base,'Results_AllConfigs.csv'),'TextType','string');
B=readtable(fullfile(base,'Baseline_PerWalk.csv'),'TextType','string');
names=["Chapter4_fixed45";"F8";"F12";"FN15";"FinalNoEdge";"Final"];
scores=nan(6,4); exact=zeros(6,1); under=exact; over=exact;
for k=1:numel(names)
    if k==1
        est=[B.Speed_mps B.Count B.StrideDuration_s B.StrideLength_m];
        ref=[B.ReferenceSpeed_mps B.ReferenceCount B.ReferenceStrideDuration_s B.ReferenceStrideLength_m];
    else
        R=T(T.config==names(k) & T.walkType=="Single pass",:);
        est=[R.estSpeed_mps R.estStepCount R.estStrideDuration_mean_s R.estStrideLength_mean_m];
        ref=[R.referenceSpeed_WAV_mps R.gtStepCount R.gtStrideDuration_s R.gtStrideLength_m];
    end
    assert(size(est,1)==142 && all(isfinite(est),'all'),'Incomplete stage results.');
    scores(k,:)=100-sqrt(mean((100*(est-ref)./ref).^2,1));
    exact(k)=sum(est(:,2)==ref(:,2));
    under(k)=sum(est(:,2)<ref(:,2));
    over(k)=sum(est(:,2)>ref(:,2));
end
stages=table(names,exact,under,over,scores(:,1),scores(:,2),scores(:,3),scores(:,4), ...
    'VariableNames',{'Stage','Exact','Under','Over','Speed_A','Count_A','Duration_A','Length_A'});
writetable(stages,fullfile(out,'Table_C2_Stages.csv'));
disp(stages);

%% Table 5.5 uses the 14 October benchmark, rather than the later rerun.
S=readtable(fullfile(base,'Sim_Trials.csv'),'TextType','string');
M=readtable(fullfile(base,'Measured_Walks.csv'),'TextType','string');
M=M(M.singlePass==1,:);
methods=["Strongest bin";"Centroid"];
fields={'errSpeed_Final_torso_peak_','errSpeed_M3_centroid'};
measuredFields={'estSpeed_Final_torso_peak_','estSpeed_M3_centroid'};
passed=zeros(2,1); pooled=passed; bias=passed; worst=passed; measured=passed; participant=nan(2,3);
for k=1:2
    e=S.(fields{k});
    pooled(k)=100-sqrt(mean(e.^2));bias(k)=mean(e);
    cells=unique([S.speed S.snr],'rows');
    a=zeros(size(cells,1),1);
    for j=1:size(cells,1)
        x=e(S.speed==cells(j,1) & S.snr==cells(j,2));
        a(j)=100-sqrt(mean(x.^2));
    end
    passed(k)=sum(a>=95);worst(k)=min(a);
    e=100*(M.(measuredFields{k})-M.gtSpeed_mps)./M.gtSpeed_mps;
    measured(k)=100-sqrt(mean(e.^2));
    for j=1:3
        x=e(M.participant=="P"+j);
        participant(k,j)=100-sqrt(mean(x.^2));
    end
end
benchmark=table(methods,passed,pooled,bias,worst,measured,participant(:,1),participant(:,2),participant(:,3), ...
    'VariableNames',{'Method','Passed','Simulated_A','Simulated_Bias_pct','Worst_A','Measured_A','P1_A','P2_A','P3_A'});
writetable(benchmark,fullfile(out,'Table_5_5_Benchmark.csv'));
disp(benchmark);

%% Table C.3 retains the Octave results used in the report.
A=readtable(fullfile(base,'ArmLegAblation_Results.csv'),'TextType','string','VariableNamingRule','preserve');
cases=unique(A.("case"),'stable');
meanBias=nan(numel(cases),1);low=meanBias;high=meanBias;edge=meanBias;
for k=1:numel(cases)
    rows=A.("case")==cases(k);
    meanBias(k)=mean(A.bias_pct(rows));low(k)=min(A.bias_pct(rows));high(k)=max(A.bias_pct(rows));
    edge(k)=mean(A.lowerEdge_rel(rows));
end
ablation=table(cases,meanBias,low,high,edge);
writetable(ablation,fullfile(out,'Table_C3_Ablation.csv'));
disp(ablation);
end
