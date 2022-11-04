clear all
close all force

addpath(genpath('../../helpers/'))

load ../data/ISC_ECG_Paris.mat


%%% figure options
FontSize = 11;
MarkerEdgeColor = 'k';
Marker = 'o';
LineColor = [1 1 1]*0.65;

%%% remove incomplete subjects

bad_subjects = [4 7 12 13];
good_subjects = setdiff(1:25, bad_subjects); %return the values of A that are not in B without repetition

SUB_ISC_CA = isc_per_subject_A(good_subjects,:,:);
SUB_ISC_CD = isc_per_subject_D(good_subjects,:,:);

clear isc*


 iscCA = tanh(nanmean(atanh(SUB_ISC_CA),3));
 iscCD = tanh(nanmean(atanh(SUB_ISC_CD),3));
     
 Nrand = size(iscCA,2);
 
 pval_pos_CA = max(mean(iscCA(:,end)<iscCA(:,1:Nrand-1),2),1/Nrand); 
 pval_neg_CA = max(mean(iscCA(:,end)>iscCA(:,1:Nrand-1),2),1/Nrand);
 pval_pos_CD = max(mean(iscCD(:,end)<iscCD(:,1:Nrand-1),2),1/Nrand); 
 pval_neg_CD = max(mean(iscCD(:,end)>iscCD(:,1:Nrand-1),2),1/Nrand);
                  
 p_fdr_C = fdr([pval_neg_CA pval_neg_CD pval_pos_CA pval_pos_CD],0.01);
     
 
 n = 0;
 
 clear y_ISC
 
 for i =1:size(iscCA,1)
   
     n = n+1;
     
     y_ISC(n) = iscCA(i,end);
     x_subject(n)    = i;
     x_attention(n)  = 1;
     x_mode(n)       = 1;
  
        n = n+1;
     
     y_ISC(n) = iscCD(i,end);
     x_subject(n)    = i;
     x_attention(n)  = 2;
     x_mode(n)       = 1;
     
     
 end
     
 fig = figure('Units','inches','Position',[10    6    10    5],'Color','White');

 
               %%%% ISC plot
         sa =     subplot(121);
              
      plot(0.85,iscCA(:,end),'ro')
          hold on
      if sum(p_fdr_C(:,3))>0
      plot(0.85,iscCA(logical(p_fdr_C(:,3)),end),'ro','MarkerFaceColor',[1 0 0])
      end
      if sum(p_fdr_C(:,1))>0
      plot(0.85,iscCA(logical(p_fdr_C(:,1)),end),'ro','MarkerFaceColor',[0 0 1])
      end
      
      plot(1.15,iscCD(:,end),'bo')
      if sum(p_fdr_C(:,4))>0
      plot(1.15,iscCD(logical(p_fdr_C(:,4)),end),'bo','MarkerFaceColor',[1 0 0])
      end
      if sum(p_fdr_C(:,2))>0
      plot(1.15,iscCD(logical(p_fdr_C(:,2)),end),'bo','MarkerFaceColor',[0 0 1])
      end    
      
      for i = 1:size(iscCD,1)
          line([0.85 1.15],[iscCA(i,end) iscCD(i,end)],'color',LineColor)
      end
      
% average relationship
line([0.85 1.15], [mean(iscCA(:,end)) mean(iscCD(:,end))],'color','k','LineWidth',2)
plot([0.85 1.15], [mean(iscCA(:,end)) mean(iscCD(:,end))],'Marker',Marker,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceColor',[0 0 0])

set(sa,'XTick',[0.85 1.15])
set(sa,'XTickLabel',{'Attentive', 'Distracted'})
set(sa,'FontSize',FontSize+2)
ylabel('ISC-HR','FontSize',FontSize+3)

sb_B_pos = get(sa,'Position'); offset_x = 0.04; offset_y = 0.06;
annotation('textbox',[sb_B_pos(1)-offset_x sb_B_pos(2)+offset_y sb_B_pos(3) sb_B_pos(4)], 'string', 'A','FontSize', FontSize+5,'LineStyle','none')


       axis([0.7 1.3 -0.06 0.15])
       

       data_stats = iscCA(:,end) - iscCD(:,end);
       
    disp('-----')
    disp('Figure 5A - ISC-HR Paris dataset attentive versus distracted')   
       
    
    [H,P,CI,STATS] = ttest(data_stats);

    
    disp(['ttest - t(' num2str(STATS.df) ') = ' num2str(STATS.tstat) ...
        'p = ' num2str(P)])
    disp(['signrank - p-value = ' num2str(signrank(data_stats))])

    
    [bf10,p] = bf.ttest(data_stats);
 
    disp(['bf10 = ' num2str(bf10)])
    disp(['bf01 = ' num2str(1/bf10)])
    disp('-----')

       

%%%%%%%
%%%%% ISC versus BEH 

load ../data/behavior.mat


sa = subplot(122);

scatter(iscCD(:,end),score_nat,'bo')
hold on
scatter(iscCA(:,end),score_at,'ro')

      scatter(iscCA(logical(p_fdr_C(:,3)),end),score_at(logical(p_fdr_C(:,3))),'ro','MarkerFaceColor',[1 0 0])

      
 [rA,pA] = corr(iscCA(:,end),score_at,'type','Spearman');
 
 
 [rD,pD] = corr(iscCD(:,end),score_nat,'type','Spearman');
 
 dfD = length(score_nat)-2;
 dfA = length(score_at)-2;
 
 

set(gca,'Box','on')
set(gca,'FontSize',FontSize+2)

sb_C_pos = get(gca,'Position');

disp('Figure 5B - Memory versus ISC-HR')

fprintf('Spearman distracted: r(%d)=%1.3f, p=%1.5f\n',dfD,rD,pD)
[bf10,r,p] = bf.corr(iscCD(:,end),score_nat);
fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)



fprintf('Spearman Attentive: r(%d)=%1.3f, p=%1.5f\n',dfA,rA,pA)
[bf10,r,p] = bf.corr(iscCA(:,end),score_at);

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)

x = [iscCA(:,end); iscCD(:,end)];
y = [score_at; score_nat];

 [rALL,pALL] = corr(x,y,'type','Spearman');
 
 dfALL = length(x)-2;
 
fprintf('Spearman All: Beh vs. ISC : r(%d)=%1.3f, p=%f\n',dfALL,rALL,pALL)
[bf10,r,p] = bf.corr(x,y);

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)
 


xlabel('ISC-HR','FontSize', FontSize+3)
ylabel('Memory performance','FontSize', FontSize+3)

axis([-.05 .15 -10 110])

set(gca,'YTick',0:25:100)

sb_B_pos = get(sa,'Position'); offset_x = 0.04; offset_y = 0.06;
annotation('textbox',[sb_B_pos(1)-offset_x sb_B_pos(2)+offset_y sb_B_pos(3) sb_B_pos(4)], 'string', 'B','FontSize', FontSize+5,'LineStyle','none')



print('../figures/Figure5.png','-dpng','-r300',fig)


disp('-----')


%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%




%%%% HR mean and STD

%%% cardio




HRM = squeeze(nanmean(HR_mean(:,:,good_subjects)));
HRV = squeeze(nanmean(HR_std(:,:,good_subjects)));


fig = figure('Units','inches','Position',[10    6    14    14],'Color','White');


sb_A = subplot(2,4,1:2); %% cardio HR
plot(0.85,HRM(1,:),Marker,'MarkerEdgeColor','r','MarkerFaceColor','r')
hold on
plot(1.15,HRM(2,:),Marker,'MarkerEdgeColor','b','MarkerFaceColor','b')
for j = 1:size(HRM,2)
line([0.85 1.15],[HRM(1,j) HRM(2,j)],'color',LineColor)
end
% average relationship
line([0.85 1.15], mean(HRM,2),'color','k','LineWidth',2)
plot([0.85 1.15], mean(HRM,2),'Marker',Marker,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceColor',[0 0 0])

set(sb_A,'XTick',[0.85 1.15])
set(sb_A,'XTickLabel',{'Attentive', 'Distracted'})
set(sb_A,'FontSize',FontSize+2)
ylabel('HR','FontSize',FontSize+3)

sb_B_pos = get(sb_A,'Position'); offset_x = 0.05; offset_y = 0.06;
annotation('textbox',[sb_B_pos(1)-offset_x sb_B_pos(2)+offset_y sb_B_pos(3) sb_B_pos(4)], 'string', 'A','FontSize', FontSize+5,'LineStyle','none')


 axis([0.7 1.3 40 120])
 
[H,P,CI,STATS] = ttest(diff(HRM));
tstat = STATS.tstat;

title('Average','FontSize',FontSize+5)

disp('Supplementary Figure 4A - Heart rate, attentive versus distracted')

data_stats = diff(HRM);

    [bf10,p] = bf.ttest(data_stats);

    disp(['t-test HRM, t(' num2str(size(HRM,2)-1) ') = ' num2str(tstat) ', p = ' num2str(P,2)])

    disp(['signrank - p-value = ' num2str(signrank(data_stats))])
    disp(['bf10 = ' num2str(bf10)])
    disp(['bf01 = ' num2str(1/bf10)])
    disp('-----')



sb_B = subplot(2,4,3:4); %% cardio HR
plot(0.85,HRV(1,:),Marker,'MarkerEdgeColor','r','MarkerFaceColor','r')
hold on
plot(1.15,HRV(2,:),Marker,'MarkerEdgeColor','b','MarkerFaceColor','b')
for j = 1:size(HRM,2)

line([0.85 1.15],[HRV(1,j) HRV(2,j)],'color',LineColor)
end
% average relationship
line([0.85 1.15], mean(HRV,2),'color','k','LineWidth',2)
plot([0.85 1.15], mean(HRV,2),'Marker',Marker,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceColor',[0 0 0])

set(sb_B,'XTick',[0.85 1.15])
set(sb_B,'XTickLabel',{'Attentive', 'Distracted'})
set(sb_B,'FontSize',FontSize+2)
ylabel('HRV','FontSize',FontSize+3)

sb_B_pos = get(sb_B,'Position'); offset_x = 0.05; offset_y = 0.06;
annotation('textbox',[sb_B_pos(1)-offset_x sb_B_pos(2)+offset_y sb_B_pos(3) sb_B_pos(4)], 'string', 'B','FontSize', FontSize+5,'LineStyle','none')

title('Average','FontSize',FontSize+5)


axis([0.7 1.3 0 14])

disp('Supplementary Figure 4B - Heart rate variability, attentive versus distracted')


[H,P,CI,STATS] = ttest(diff(HRV));
tstat = STATS.tstat;

disp(['t-test HRV, t(' num2str(size(HRV,2)-1) ') = ' num2str(tstat) ', p = ' num2str(P,2)])

data_stats = diff(HRV);

    [bf10,p] = bf.ttest(data_stats);

    disp(['signrank - p-value = ' num2str(signrank(data_stats))])
    disp(['bf10 = ' num2str(bf10)])
    disp(['bf01 = ' num2str(1/bf10)])
    disp('-----')



%%%% correlations cardio

xlims = [0 0.15];

s_A = subplot(245);
scatter(iscCA(:,end),HRM(1,:),200,'r.')
[r,p,df,lh,ph] = ...
    line_function_confidence(iscCA(:,end),...
    HRM(1,:)',0.05,[xlims(1)+0.001 xlims(2)-0.001]);
sb_A_pos = get(s_A,'Position');
ah = annotation_paper([sb_A_pos(1)+0.07 sb_A_pos(2)+sb_A_pos(4)-0.06 sb_A_pos(3)/3+0.025 0.05],r,p,df,FontSize+3);
set(ah,'BackgroundColor','White')

annotation('textbox',[sb_A_pos(1)-offset_x sb_A_pos(2)+offset_y sb_A_pos(3) sb_A_pos(4)], 'string', 'C','FontSize', FontSize+5,'LineStyle','none')


ylabel('HR','FontSize',  FontSize+3)
xlabel('ISC-HR','FontSize', FontSize+3)
set(gca,'FontSize',FontSize+2)

title('Attentive','FontSize',FontSize+5)
set(s_A,'Box','on')
xlim(xlims)
ylim([40 120])


disp('Supplementary Figure 4C - Correlation Heart rate versus ISC-HR (attentive)')

fprintf('Pearson: r(%d)=%1.3f, p=%1.5f\n',df,r,p)
[bf10,r,p] = bf.corr(iscCA(:,end),HRM(1,:)');

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)
disp('-')

s_C = subplot(247);
scatter(iscCA(:,end),HRV(1,:),200,'r.')
[r,p,df,lh,ph] = ...
    line_function_confidence(iscCA(:,end),...
    HRV(1,:)',0.05,[xlims(1)+0.001 xlims(2)-0.001]);
sb_A_pos = get(s_C,'Position');
ah = annotation_paper([sb_A_pos(1)+0.07 sb_A_pos(2)+sb_A_pos(4)-0.06 sb_A_pos(3)/3+0.025 0.05],r,p,df,FontSize+3);
set(ah,'BackgroundColor','White')

annotation('textbox',[sb_A_pos(1)-offset_x sb_A_pos(2)+offset_y sb_A_pos(3) sb_A_pos(4)], 'string', 'D','FontSize', FontSize+5,'LineStyle','none')


ylabel('HRV','FontSize',  FontSize+3)
xlabel('ISC-HR','FontSize', FontSize+3)
set(gca,'FontSize',FontSize+2)
set(s_C,'Box','on')
xlim(xlims)
ylim([0 14])

title('Attentive','FontSize',FontSize+5)


xlims = [-0.04 0.06];

disp('Supplementary Figure 4D - Correlation Heart rate variability versus ISC-HR (attentive)')

fprintf('Pearson: r(%d)=%1.3f, p=%1.5f\n',df,r,p)
[bf10,r,p] = bf.corr(iscCA(:,end),HRV(1,:)');

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)
disp('-')


s_B = subplot(246);
scatter(iscCD(:,end),HRM(2,:),200,'b.')
[r,p,df,lh,ph] = ...
    line_function_confidence(iscCD(:,end),...
    HRM(2,:)',0.05,[xlims(1)+0.001 xlims(2)-0.001]);
sb_A_pos = get(s_B,'Position');
ah = annotation_paper([sb_A_pos(1)+0.07 sb_A_pos(2)+sb_A_pos(4)-0.06 sb_A_pos(3)/3+0.025 0.05],r,p,df,FontSize+3);
set(ah,'BackgroundColor','White')


set(gca,'YTickLabel',[])


xlabel('ISC-HR','FontSize', FontSize+3)
set(gca,'FontSize',FontSize+2)
set(s_B,'Box','on')
title('Distracted','FontSize',FontSize+5)
xlim(xlims)
ylim([40 120])


disp('Supplementary Figure 4C - Correlation Heart rate versus ISC-HR (distracted)')

fprintf('Pearson: r(%d)=%1.3f, p=%1.5f\n',df,r,p)
[bf10,r,p] = bf.corr(iscCD(:,end),HRM(2,:)');

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)
disp('-')

s_D = subplot(248);
scatter(iscCD(:,end),HRV(2,:),200,'b.')
[r,p,df,lh,ph] = ...
    line_function_confidence(iscCD(:,end),...
    HRV(2,:)',0.05,[xlims(1)+0.001 xlims(2)-0.001]);
sb_A_pos = get(s_D,'Position');
ah = annotation_paper([sb_A_pos(1)+0.07 sb_A_pos(2)+sb_A_pos(4)-0.06 sb_A_pos(3)/3+0.025 0.05],r,p,df,FontSize+3);
set(ah,'BackgroundColor','White')



xlabel('ISC-HR','FontSize', FontSize+3)
set(gca,'FontSize',FontSize+2)
set(gca,'YTickLabel',[])

set(s_D,'Box','on')
xlim(xlims)
%print(sprintf('../figures/ISC_HRV_by_condition_cardio.png'),'-dpng','-r300',fig)
ylim([0 14])


title('Distracted','FontSize',FontSize+5)


disp('Supplementary Figure 4D - Correlation Heart rate variability versus ISC-HR (distracted)')

fprintf('Pearson: r(%d)=%1.3f, p=%1.5f\n',df,r,p)
[bf10,r,p] = bf.corr(iscCD(:,end),HRV(2,:)');

fprintf('b10=%1.2f, bf01=%1.2f\n',bf10,1/bf10)
disp('-')

print('../figures/Sup_Figure4.png','-dpng','-r300',fig)

