clearvars
close all
addpath(genpath('../../helpers'))
addpath ./helpers/

%%first load the raw ECG data
options.directory_data = '../data/raw/';
options.directory_metadata = '../data/metadata/';
options.directory_results_csv = '../results/frequency/';

%% figure options
FontSize = 10;
MarkerEdgeColor = 'k';
Marker = 'o';
ylim_individual = [-0.05 0.3];
ylim_average    = [-0.05 0.3];
xlim_offset = 0.5;
LineColor = [1 1 1]*0.65;

color_attend = 'Red';
color_distract = 'Blue';

%% HR paramaters
hr_int_method = 'pchip';
fs = 2048;
fs_new = 128;

% number of circular shifts for permutation stats
Nrand=10000;
seed = 123456;
rng(seed,'twister')

%% filter paramaters
N = 2^16; %points for FFT
filter_order = 5;
center_frequencies = 2.^(-5:0.1:2);
Nfrequency = length(center_frequencies);

%% create and plot filters magnitude response
for iFrequency = 1:Nfrequency
    % pick your preferred band-pass filter
    cutoff_frequencies(iFrequency,:) = [0.9 1.1]*center_frequencies(iFrequency);
    
    [Z,P,k]=butter(filter_order,cutoff_frequencies(iFrequency,:)/(fs/2),'bandpass' );    
    [sos{iFrequency,1}] = zp2sos(Z,P,k);
    
end

%% load the data, preprocess and compute HR

% first load the ECG data
[data_ecg,time_axis,metadata] = loadNYdata(options);

Nsubjects = size(data_ecg,1);

for iSubject = 1:Nsubjects
    
    fprintft('Processing data for participant %d/%d (%d)\n',iSubject,Nsubjects,metadata{iSubject}.participant_no) %
    %preprocess the data (i.e. filtering etc.)
    tmp = preprocessECGdata(data_ecg{iSubject},fs);

    %compute Heart rate from the data
    data_hr{iSubject,1} = computeHRdata(tmp,fs,hr_int_method);
    
end
clear data_ecg

for iFrequency = 1:Nfrequency
    
    % filter data
    fprintf('Bandpass with center frequency: %1.5fHz (%d/%d)\n',center_frequencies(iFrequency),iFrequency,Nfrequency)
    for iSubject = 1:Nsubjects
        data_hr_bp{iSubject,1} = sosfilt(sos{iFrequency,1},data_hr{iSubject,1}-data_hr{iSubject,1}(1),1);
        
        if any(isnan(data_hr_bp{iSubject,1}))
            figure
            plot(data_hr_bp{iSubject,1})
            title(sprintf('Subject=%d, freq_no=%d',iSubject,iFrequency))
            drawnow
        end
        
        %segment the data
        data_hr_bp{iSubject,1} = segmentData(data_hr_bp{iSubject,1},metadata{iSubject,1},'exg','stimuli');
        fprintf('.')
    end
    
    fprintf('done\n')
    
    %aggregate all heart rate data and truncate ends if needed
    data_hr_bp = aggregateData(data_hr_bp,fs)';

    %resample the data (from fs -> fs_new)
    data_hr_freq{iFrequency} = resampleEXGdata(data_hr_bp,fs,fs_new);
end
clear data_hr_bp data_hr

%% start ISC of Heart rate computations
Nstim  = size(data_hr_freq{1},2);
Nconditions     = size(data_hr_freq{1},1);

for iFrequency = 1:Nfrequency
fprintf('Center frequency: %1.5fHz',center_frequencies(iFrequency))
    
    for iStim=1:Nstim
        for iCondition=1:Nconditions
            fprintf('.')
           
            %get the heart rate data
            HR = data_hr_freq{iFrequency}{iCondition,iStim};
            HR_attend = data_hr_freq{iFrequency}{1,iStim};
            
            [T,Nsubjects] = size(HR);
            [T_attend,Nsubjects_attend] = size(HR_attend);
            
            %compute ISC
            ISC = corr(HR,HR_attend); ISC(1:Nsubjects+1:Nsubjects^2)=NaN; %set diaginal to NaN
            ISC_per_subject(:,iStim,iCondition,iFrequency) = tanh(nanmean(atanh(ISC))); %average across subjects

            HRV_per_subject(:,iStim,iCondition,iFrequency) = sqrt(mean((HR-mean(HR(:))).^2,1));
        end
    end
    fprintf('done\n')
end

ISC_per_frequency = squeeze(tanh(mean(atanh(ISC_per_subject),2)));
HRV_per_frequency = squeeze(mean(HRV_per_subject,2));

ISC_per_frequency_subject_average = squeeze(tanh(mean(mean(atanh(ISC_per_subject),2),1)));
HRV_per_frequency_subject_average = squeeze(mean(mean(HRV_per_subject,2),1));

HRV_per_frequency_subject_stderr = squeeze(std(HRV_per_frequency,1))./sqrt(Nsubjects);
ISC_per_frequency_subject_stderr = squeeze(atanh(std(tanh(ISC_per_frequency),1)))./sqrt(Nsubjects);

%compute uncorrected paired t-test between the two conditions for each frequency bin
for iFrequency = 1:Nfrequency
   delta_ISC(iFrequency,:) = squeeze(ISC_per_frequency(:,1,iFrequency)-ISC_per_frequency(:,2,iFrequency));
   [~,pval_isc_attend_distract(iFrequency)] = ttest(delta_ISC(iFrequency,:));   
   
   delta_HRV(iFrequency,:) = squeeze(HRV_per_frequency(:,1,iFrequency)-HRV_per_frequency(:,2,iFrequency));
   [~,pval_hrv_attend_distract(iFrequency)] = ttest(delta_HRV(iFrequency,:));
end

% %% write results to csv files
% for iFrequency = 1:Nfrequency
%     for iCondition=1:Nconditions
%         filename_isc_frequency = sprintf('ISC_HR_frequency_%s_allstim_condition=%d.csv',strrep(sprintf('%1.4f',center_frequencies(iFrequency)),'.','_'),iCondition);
%         
%         if ~exist([options.directory_results_csv filename_isc_frequency],'file')
%             writematrix(ISC_per_subject(:,:,iCondition,iFrequency),[options.directory_results_csv filename_isc_frequency])
%         end
%         
%         filename_hrv_frequency = sprintf('HRV_frequency_%s_allstim_condition=%d.csv',strrep(sprintf('%1.4f',center_frequencies(iFrequency)),'.','_'),iCondition);
%         
%         if ~exist([options.directory_results_csv filename_hrv_frequency],'file')
%             writematrix(HRV_per_subject(:,:,iCondition,iFrequency),[options.directory_results_csv filename_hrv_frequency])
%         end
%     end
% end

%% do cluster permutation test

% cluster stats parameters
cfg.thr1 = 0.05;
cfg.thr2 = 0.05;
cfg.method = 'tstat'; % 'wcr' = ranksum, 'tstat' = tstat / ranksum is much slower 
cfg.n_perm = 10000;


data_stats = permute(HRV_per_frequency, [3 1 2]);

[clusters_hrv_attend_distract, p_value] = clust_1d_paired(data_stats, cfg);
significant_cluster_hrv_attend_distract = clusters_hrv_attend_distract>0;

clusters = unique(clusters_hrv_attend_distract);

disp('Clusters HRV')
disp('------------')

for j = 1:(length(clusters)-1)

    index = find(clusters_hrv_attend_distract==clusters(j+1));
    
    f_min = center_frequencies(min(index));
    f_max = center_frequencies(max(index));
    
    disp(['Cluster #' num2str(j) ' for ' num2str(f_min) ' Hz to ' num2str(f_max) ' Hz. p-value < ' num2str(p_value(j))])
    
end

disp('=')



data_stats = permute(ISC_per_frequency, [3 1 2]);

[clusters_isc_attend_distract, p_value] = clust_1d_paired(data_stats, cfg);
significant_cluster_isc_attend_distract = clusters_isc_attend_distract>0;

disp('Clusters ISC')
disp('------------')

clusters = unique(clusters_isc_attend_distract);

for j = 1:(length(clusters)-1)

    index = find(clusters_isc_attend_distract==clusters(j+1));
    
    f_min = center_frequencies(min(index));
    f_max = center_frequencies(max(index));
    
    disp(['Cluster #' num2str(j) ' for ' num2str(f_min) ' Hz to ' num2str(f_max) ' Hz. p-value < ' num2str(p_value(j))])
    
end

%% plot average across all stories and subjects with condidence intervals (cluster statistics)
fig = figure('units','inches','position',[10    6    7    3.5],'Color','White');
sb_A = subplot(1,2,1);

% panel A
[s_attend,ph_attend] = boundedline(center_frequencies, HRV_per_frequency_subject_average(1,:), HRV_per_frequency_subject_stderr(1,:),'-k.','transparency',0.05);
set(s_attend,'LineWidth',1), set(s_attend,'Marker','none'), hold on, set(s_attend,'Color',color_attend),set(ph_attend,'FaceColor',color_attend),set(ph_attend,'FaceAlpha',0.2)

[s_distracted,ph_distracted] = boundedline(center_frequencies, HRV_per_frequency_subject_average(2,:), HRV_per_frequency_subject_stderr(2,:),'-r.','transparency',0.05);
set(s_distracted,'LineWidth',1), set(s_distracted,'Marker','none'), set(s_distracted,'Color',color_distract),set(ph_distracted,'FaceColor',color_distract),set(ph_distracted,'FaceAlpha',0.2)

set(gca,'Xscale','log')
xlabel('Frequency [Hz]')
ylabel('HRV [BPM]')
set(gca,'FontSize',FontSize)
xlim([min(center_frequencies) max(center_frequencies)])
set(sb_A,'Xtick',10.^(-4:1:0))
set(sb_A,'Position',get(sb_A,'Position')+[-0.05 0.03 0.07 0])
plot_significance_areas(sb_A,significant_cluster_hrv_attend_distract,[],cutoff_frequencies,get(sb_A,'Ylim'))

sb_A_pos = get(sb_A,'Position'); offset_x = 0.07; offset_y = 0.06;
annotation('textbox',[sb_A_pos(1)-offset_x sb_A_pos(2)+offset_y sb_A_pos(3) sb_A_pos(4)], 'string', 'A','FontSize', FontSize+3,'LineStyle','none')
uistack(s_attend,'top')
uistack(s_distracted,'top')
set(gca,'Box','on')

% panel B
sb_B = subplot(1,2,2);
[s_attend,ph_attend] = boundedline(center_frequencies, ISC_per_frequency_subject_average(1,:), ISC_per_frequency_subject_stderr(1,:),'-k.','transparency',0.05);
set(s_attend,'LineWidth',1), set(s_attend,'Marker','none'), hold on, set(s_attend,'Color',color_attend),set(ph_attend,'FaceColor',color_attend),set(ph_attend,'FaceAlpha',0.2)

[s_distracted,ph_distracted] = boundedline(center_frequencies, ISC_per_frequency_subject_average(2,:), ISC_per_frequency_subject_stderr(2,:),'-r.','transparency',0.05);
set(s_distracted,'LineWidth',1), set(s_distracted,'Marker','none'), set(s_distracted,'Color',color_distract),set(ph_distracted,'FaceColor',color_distract),set(ph_distracted,'FaceAlpha',0.2)

set(gca,'Xscale','log')
xlabel('Frequency [Hz]')
ylabel('ISC-HR')
set(gca,'FontSize',FontSize)
xlim([min(center_frequencies) max(center_frequencies)])
ylim([-0.015 0.06])
set(gca,'Xtick',10.^(-4:1:0))
set(gca,'Ytick',(-1:5)/1e2)
set(sb_B,'Position',get(sb_B,'Position')+[0.01 0.03 0.07 0])
plot_significance_areas(sb_B,significant_cluster_isc_attend_distract,[],cutoff_frequencies,get(sb_B,'Ylim'))
legend([s_attend s_distracted],{'Attentive','Distracted'},'FontSize',FontSize)

sb_B_pos = get(sb_B,'Position'); offset_x = 0.07; offset_y = 0.06;
annotation('textbox',[sb_B_pos(1)-offset_x sb_B_pos(2)+offset_y sb_B_pos(3) sb_B_pos(4)], 'string', 'B','FontSize', FontSize+3,'LineStyle','none')
uistack(s_attend,'top')
uistack(s_distracted,'top')
set(gca,'Box','on')

print(sprintf('../figures/Figure4.png'),'-dpng','-r300',fig)



%% plot average across subjects

ISC_per_story = tanh(squeeze(mean(atanh(ISC_per_subject),1)));

colors = {'Green','Orange','Purple','Brown','Gold'};

fig=figure;
for iStim = 1:Nstim
    s_attend(iStim) = semilogx(center_frequencies,squeeze(ISC_per_story(iStim,1,:)),'Color',rgb(colors{iStim}),'LineWidth',2,'LineStyle','-'); hold on
    s_distracted(iStim) = semilogx(center_frequencies,squeeze(ISC_per_story(iStim,2,:)),'Color',rgb(colors{iStim}),'LineWidth',2,'LineStyle','--');
end
legend(s_attend,cellfun(@(s) sprintf('Story %d',s),num2cell(1:5),'UniformOutput',false),'FontSize',FontSize)
xlabel('Center frequency [Hz]')
ylabel('ISC-HR')
set(gca,'FontSize',FontSize)
xlim([min(center_frequencies) max(center_frequencies)])
set(gca,'Xtick',10.^(-4:1:0))

%print(sprintf('../figures/ISC_frequency_resolved_by_story.png'),'-dpng','-r300',fig)

