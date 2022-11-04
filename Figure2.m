clc

clearvars
close all
addpath(genpath('../../helpers/'))

installBayesFactor

%options.directory_segmented_data_csv = '../data/segmented/csv/';
%options.directory_results_csv = '../results/';

%%% set the number of permutations for testing
Npermutations=10000;

%load data , sampling frequency 25 hz 
%%% 27 subjects, 16 recording segments 
load ../data/BHX_fs25.mat

%set the seed
seed = 123456;
rng(seed,'twister')

%set the corr method
cmet = 'Pearson';%'Pearson' 

%% plot settings
FontSize = 10;

%% down sample heart rate data
Nsegments = length(HR_interp);
fs_new = 4; %% new sampling frequency == 4 hz 

for iSegment=1:Nsegments
    valid_index = find(~isnan(sum(HR_interp{iSegment},2))); % remove NaN (buggy if NaN not consecutive)
    HR_interp{iSegment} = resample(HR_interp{iSegment}(valid_index,:)-HR_interp{iSegment}(valid_index(1),:),fs_new,fs)+HR_interp{iSegment}(valid_index(1),:); % slow signal anyway
end

% %% write out segmented/preprocessed data in to CSV files
% for iSegment=1:Nsegments
%     
%     filename_segmented_csv = sprintf('HR_instantaneous_segment=%d.csv',iSegment);
%     if ~exist([options.directory_segmented_data_csv filename_segmented_csv],'file')
%         writematrix(HR_interp{iSegment},[options.directory_segmented_data_csv filename_segmented_csv])
%     end
% end

% %% show me the raw data
% fig = figure('units','normalized','Position',[0.25    0.16    0.5    0.67],'Color','White');
% for iSegment=1:Nsegments
%     subplot(4,4,iSegment)
%     HR = HR_interp{iSegment};
%     HR = HR - mean(HR);
%     [T,N] = size(HR);
%     imagesc((1:T)/fs_new,1:N,HR'); caxis([-20 20])
%     
%     if rem(iSegment-1,4)~=0, set(gca,'ytick',[]); else, ylabel('Subject no.','FontSize',FontSize); end
%     if iSegment<=12,         set(gca,'xtick',[]); else, xlabel('Time [s]','FontSize',FontSize); end
%     
%     title(sprintf('Segment %d',iSegment),'FontSize',FontSize)
%     set(gca,'Ydir','normal')
% end
% 
% title('Heart rate','FontSize',FontSize+5)
% 
% %% show all correlations
% fig = figure('units','normalized','Position',[0.25    0.16    0.42    0.7],'Color','White');
% for iSegment=1:Nsegments
%     h = subplot(4,4,iSegment);
%     HR = HR_interp{iSegment};
%     [T,N] = size(HR);
% 
%         ISC = corr(HR,'type',cmet); ISC(1:N+1:N^2)=NaN;
%     imagesc(ISC); %caxis([-.3 0.3])
%     title(sprintf('Segment %d',iSegment),'FontSize',FontSize)
%     if rem(iSegment-1,4)~=0, set(h,'ytick',[]); else, ylabel('Subject no.','FontSize',FontSize); end
%     if iSegment<=12,         set(h,'xtick',[]); else, xlabel('Subject no.','FontSize',FontSize); end
%     set(gca,'Ydir','normal')
%     set(gca,'DataAspectRatio',[1 1 1])
% end
% pos = get(h,'Position'); hbar = colorbar; set(h,'Position',pos)
% title(hbar,'r')
% 
% title('ISC of heart rate','FontSize',FontSize+5)

%% A: ISC per subject and shuffle stats

for iSegment=1:Nsegments
    disp(['segment ' num2str(iSegment) ' of ' num2str(Nsegments)])
    HR = zscore(HR_interp{iSegment});
    [T,Nsubjects] = size(HR);

       ISC = atanh(corr(HR,'type', cmet)); ISC(1:Nsubjects+1:Nsubjects^2)=NaN;
       ISC_per_subject(:,iSegment) = tanh(nanmean(ISC));
    
    HR_shifted = HR;
    
    %compute ISC based on random circular shifts
    for iRandShift=Npermutations:-1:1
        % circular shuffle, no tricks.
        for iSubject=1:Nsubjects, HR_shifted(:,iSubject) = circshift(HR(:,iSubject),round(T*rand(1))); end
        
        %compute ISC
        
        ISC = atanh(corr(HR_shifted,HR_shifted,'type',cmet)); ISC(1:Nsubjects+1:Nsubjects^2)=NaN; %set diaginal to NaN
        ISC_per_subject_shifted(:,iRandShift,iSegment) = tanh(nanmean(ISC)); %average across subjects
    end
    
    for iSubject=1:Nsubjects
        pval_pos(iSubject,iSegment) = max(mean(ISC_per_subject(iSubject,iSegment)<ISC_per_subject_shifted(iSubject,1:Npermutations,iSegment)),1/Npermutations);
        pval_neg(iSubject,iSegment) = max(mean(ISC_per_subject(iSubject,iSegment)>ISC_per_subject_shifted(iSubject,1:Npermutations,iSegment)),1/Npermutations);
        
        %if ISC values are NaN discard the stats
        if isnan(ISC_per_subject(iSubject,iSegment))
            pval_pos(iSubject,iSegment) = NaN;
            pval_neg(iSubject,iSegment) = NaN;
        end
    end
end

%determine which subject had significant ISC values for each segment using FDR
h = fdr([pval_pos pval_neg],0.01);
sig_pos = h(:,1:Nsegments)>0;
sig_neg = h(:,Nsegments+1:end)>0;

% %% write out csv files with results (ISC, HR and HRV)
% for iSegment=1:Nsegments
%     filename_isc = sprintf('ISC_HR_segment=%d.csv',iSegment);
%     filename_isc_permuted = sprintf('ISC_HR_permuted_segment=%d.csv',iSegment);
%     
%     if ~exist([options.directory_results_csv filename_isc],'file')
%         writematrix(ISC_per_subject(:,iSegment),[options.directory_results_csv filename_isc])
%     end
%     
%     if ~exist([options.directory_results_csv filename_isc_permuted],'file')
%         writematrix(ISC_per_subject_shifted(:,:,iSegment),[options.directory_results_csv filename_isc_permuted])
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the results %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 - A % %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('units','inches','position',[10    6    7    3.5],'Color','White');
sb_A = axes('Position',[0.1    0.14    0.22    0.82]);

for iSegment=1:Nsegments
    plot(iSegment,ISC_per_subject(:,iSegment),'.','MarkerEdgeColor',[1 1 1]*0.65,'MarkerSize',10); hold on
    if sum(sig_pos(:,iSegment)), plot(sb_A,iSegment,ISC_per_subject(sig_pos(:,iSegment),iSegment),'o','MarkerFaceColor','Black','MarkerEdgeColor','Black','MarkerSize',4); end
    if sum(sig_neg(:,iSegment)), plot(sb_A,iSegment,ISC_per_subject(sig_neg(:,iSegment),iSegment),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',4); end
    xlim([0 length(HR_interp)+1])
    xlabel('segment number')
    ylabel('ISC-HR')
end

%make sure the plot is still pretty
xlim(sb_A,[0 length(HR_interp)+1])
xlabel(sb_A,'Segment no.','FontSize',FontSize)
ylabel(sb_A,'ISC-HR','FontSize',FontSize)
set(sb_A,'FontSize',FontSize)
sb_A_pos = get(sb_A,'Position'); offset_x = 0.09; offset_y = 0.06;
annotation('textbox',[sb_A_pos(1)-offset_x sb_A_pos(2)+offset_y sb_A_pos(3) sb_A_pos(4)], 'string', 'A','FontSize', FontSize+3,'LineStyle','none')



%% B: show me the summary averaged over segments

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 - B % %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


sb_B = axes('Position',[0.44    0.14    0.22    0.82]);
plot(sb_B,[0 Nsubjects+1],[0 0],'k'); hold on
ISC_per_subject_mean = tanh(mean(atanh(ISC_per_subject),2));
ISC_per_subject_shifted_mean = tanh(mean(atanh(ISC_per_subject_shifted),3));

for iSubject=1:Nsubjects % 
    pval_pos(iSubject) = max(mean(ISC_per_subject_mean(iSubject)<ISC_per_subject_shifted_mean(iSubject,:)),1/Npermutations);
end

sig_pos = fdr(pval_pos,0.05)>0;
[m,iSegment] = sort(ISC_per_subject_mean); % uggly indexing coming up!!!
plot(m,'o','MarkerEdgeColor',[1 1 1]*0.65); hold on
n=find(sig_pos(iSegment)); if ~isempty(n), plot(sb_B,n,m(n),'o','MarkerFaceColor','Black','MarkerEdgeColor','Black','MarkerSize',5); end
xlim(sb_B,[0 Nsubjects+1])
ylabel(sb_B,'mean ISC-HR','FontSize',FontSize)
xlabel(sb_B,'Subject no.','FontSize',FontSize)
set(sb_B,'FontSize',FontSize)
sb_B_pos = get(sb_B,'Position'); offset_x = 0.09; offset_y = 0.06;
annotation('textbox',[sb_B_pos(1)-offset_x sb_B_pos(2)+offset_y sb_B_pos(3) sb_B_pos(4)], 'string', 'B','FontSize', FontSize+3,'LineStyle','none')

ylim([-0.015 .06])


%%%% display results 

[bf10,p] = bf.ttest(ISC_per_subject_mean);

display('figure 2 - panel B')

display(['t-test, p=' num2str(p)])

display(['Bayes, bf10=' num2str(bf10)])

p_s = signrank(ISC_per_subject_mean);

display(['signrank, p=' num2str(p_s)])


%% C: permutation of history segments

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figure 2 - C % %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear ISC

% parameters
Nstory = size(HR_interp,2);
Nsub = size(HR_interp{1},2);


clear L
%%%% equate lenght across history segments
for iStory = 1:Nstory
    L(iStory) = size(HR_interp{iStory},1);
end
% keep the mininum length
L = min(L);

for iStory = 1:Nstory
    HR_interp_fix{iStory} = zscore(HR_interp{iStory}(1:L,:));
end


% no permutation results

for iStory = 1:Nstory
        
    C = atanh(corr(HR_interp_fix{iStory},'type',cmet));
    C(1:(Nsub+1):Nsub^2) = nan;
    ISC(:,iStory) = tanh(nanmean(C));

end

% permutations over different history segments

for iRep =1:Npermutations
    aux = [];
    
%    disp(['Rep '  num2str(iRep) ' of ' num2str(Npermutations)])
    
    %assign random segments to each subject
    for iSubject = 1:Nsub
        iStory = randi(Nstory);
        aux = cat(2,aux,HR_interp_fix{iStory}(:,iSubject));
    end
    
    C = atanh(corr(aux,'type',cmet));
    C(1:(Nsub+1):Nsub^2) = nan;
    
    ISC_permut(:,iRep) = tanh(nanmean(C));
    
end


%%%% group results
mean_ISC = tanh(nanmean(atan(ISC),2));
mean_ISC_permut = tanh(nanmean(atanh(ISC_permut),2));

%% plot the results
sb_C = axes('Position',[0.77    0.14    0.22    0.82]);

plot(1,mean_ISC,'o','MarkerEdgeColor','k'), hold on
plot(2,mean_ISC_permut,'o','MarkerEdgeColor','k')


for j = 1:size(ISC,1)
    line([1 2],[mean_ISC(j) mean_ISC_permut(j)],'Color',[1 1 1]*0.65)
end

xlim([0.75 2.25])
ylim([-0.015 0.065])


ylabel('Mean ISC-HR','FontSize', FontSize)
set(gca,'Xtick',[1 2])
set(gca,'XTickLabel',{'Original', 'Permuted'})
set(gca,'FontSize',FontSize)

fprintf('--\n')


disp('Figure 2 - Panel C - Original versus permuted data')


p_group = signrank(mean_ISC-mean_ISC_permut);


[bf10,p] = bf.ttest(mean_ISC-mean_ISC_permut);


  [H,P,CI,STATS] = ttest(mean_ISC-mean_ISC_permut);

display(['paired t-test, t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ', p=' num2str(P)])

display(['Bayes, bf10=' num2str(bf10)])

display(['Bayes, bf01=' num2str(1/bf10)])

disp(['paired signrank, p-value: ' num2str(p_group)])

  [H,P,CI,STATS] = ttest(mean_ISC-mean_ISC_permut);


sb_C_pos = get(sb_C,'Position'); offset_x = 0.09; offset_y = 0.06;
annotation('textbox',[sb_C_pos(1)-offset_x sb_C_pos(2)+offset_y sb_C_pos(3) sb_C_pos(4)], 'string', 'C','FontSize', FontSize+3,'LineStyle','none')


print('../figures/Figure2.png','-dpng','-r300',fig)


%% suplementary 

axes_width_adjusment = 0.07;

factor_subject = repmat((1:Nsubjects),[Nsegments,1]);
factor_stimuli = repmat((1:Nsegments)',[1,Nsubjects]);

interactions = eye(2); %only main fixed effects
[pval, tbl, stats] = anovan(reshape(HR_mean(:),[],1),{factor_subject(:) factor_stimuli(:)}, ...
    'model', interactions, 'random',1,  'varnames', {'Subject', 'Segments'},'display','off');

fprintf('-----\n')

fprintf('Supplementary Figure 1\n')
fprintf('\n')

fprintf('ANOVA: HR ~ History segments + subjects (random)\n')
displayAnova(tbl)

interactions = eye(2); %only main fixed effects
[pval, tbl, stats] = anovan(reshape(HR_std(:),[],1),{factor_subject(:) factor_stimuli(:)}, ...
    'model', interactions, 'random',1,  'varnames', {'Subject', 'Segments'},'display','off');

fprintf('ANOVA: HRV ~ History segments + subjects (random)\n')
displayAnova(tbl)

%%
fig = figure('units','inches','position',[10    6    7    6.5],'Color','White');

% panel A
sb_A = subplot(221);

h = notBoxPlot(HR_mean');
arrayfun(@(s) set(s.data,'Marker','.'),h)
% title(['ANOVA F(X,X) = X,  p = '])
xlabel('History segments','FontSize',FontSize)
ylabel('HR','FontSize',FontSize)
set(gca,'FontSize',FontSize)
set(gca,'Box','on')
set(sb_A,'Position',get(sb_A,'Position') + [-0.05 0 axes_width_adjusment 0])

sb_A_pos = get(sb_A,'Position'); offset_x = 0.08; offset_y = 0.06;
annotation('textbox',[sb_A_pos(1)-offset_x sb_A_pos(2)+offset_y sb_A_pos(3) sb_A_pos(4)], 'string', 'A','FontSize', FontSize+3,'LineStyle','none')

% panel B
sb_B = subplot(222);
h = notBoxPlot(HR_std');
arrayfun(@(s) set(s.data,'Marker','.'),h)
xlabel('History segments','FontSize',FontSize)
ylabel('HRV','FontSize',FontSize)
set(gca,'FontSize',FontSize)
set(gca,'Box','on')
set(sb_B,'Position',get(sb_B,'Position') + [0.01 0 axes_width_adjusment 0])

sb_B_pos = get(sb_B,'Position'); offset_x = 0.08; offset_y = 0.06;
annotation('textbox',[sb_B_pos(1)-offset_x sb_B_pos(2)+offset_y sb_B_pos(3) sb_B_pos(4)], 'string', 'B','FontSize', FontSize+3,'LineStyle','none')

% panel C
xlims = ([-0.005 0.059]);
sb_C = subplot(223);

scatter(ISC_per_subject_mean,mean(HR_mean,1)',100,'k.')
xlabel('ISC-HR','FontSize', FontSize)
ylabel('HR','FontSize', FontSize)
[r, p, df, lh, ph] = line_function_confidence(ISC_per_subject_mean,mean(HR_mean,1)', 0.05,[xlims(1)+0.0001 xlims(2)-0.0001]);
set(gca,'Box','on')
set(gca,'FontSize',FontSize)
xlim(xlims)
set(sb_C,'Position',get(sb_C,'Position') + [-0.05 0 axes_width_adjusment 0])
sb_C_pos = get(sb_C,'Position'); offset_x = 0.08; offset_y = 0.06;
annotation('textbox',[sb_C_pos(1)-offset_x sb_C_pos(2)+offset_y sb_C_pos(3) sb_C_pos(4)], 'string', 'C','FontSize', FontSize+3,'LineStyle','none')

sb_C_pos = get(sb_C,'Position');
ah = annotation_paper([sb_C_pos(1)+0.01 sb_C_pos(2)+sb_C_pos(4)-0.08 sb_C_pos(3)/3 0.07],r,p,df,FontSize);
set(ah,'BackgroundColor','White')

fprintf('--\n')

fprintf('Correlation ISC vs. HR\n')

fprintf('r(%d)=%1.2f, p=%1.2f\n',df,r,p)

[bf10,r,p] = bf.corr(ISC_per_subject_mean,mean(HR_mean,1)');

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)



% panel D
sb_D = subplot(224);

scatter(ISC_per_subject_mean,mean(HR_std,1)',100,'k.')
xlabel('ISC-HR','FontSize', FontSize)
ylabel('HRV','FontSize', FontSize)
[r, p, df, lh, ph] = line_function_confidence(ISC_per_subject_mean,mean(HR_std,1)', 0.05,[xlims(1)+0.0001 xlims(2)-0.0001]);
set(gca,'Box','on')
set(gca,'FontSize',FontSize)
xlim(xlims)
set(sb_D,'Position',get(sb_D,'Position') + [0.01 0 axes_width_adjusment 0])

sb_D_pos = get(sb_D,'Position'); offset_x = 0.08; offset_y = 0.06;
annotation('textbox',[sb_D_pos(1)-offset_x sb_D_pos(2)+offset_y sb_D_pos(3) sb_D_pos(4)], 'string', 'D','FontSize', FontSize+3,'LineStyle','none')

sb_D_pos = get(sb_D,'Position');
ah = annotation_paper([sb_D_pos(1)+0.01 sb_D_pos(2)+sb_D_pos(4)-0.08 sb_D_pos(3)/3 0.07],r,p,df,FontSize);
set(ah,'BackgroundColor','White')

fprintf('--\n')

fprintf('Correlation ISC vs. HRV\n')

fprintf('r(%d)=%1.2f, p=%1.2f\n',df,r,p)

[bf10,r,p] = bf.corr(ISC_per_subject_mean,mean(HR_std,1)');

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)


print('../figures/Sup_Figure1.png','-dpng','-r300',fig)


