clc;close all;clear

%% EXCLUDE SUBJECT 9 THEIR DATA WAS OVERWRITTEN WITH THE DATA OF PARTICIPANT 6
%% Add to path Data folder (here extracted ECG folder)
addpath('./ECG_Seteb_extracted/');

%% Read files

% Read name of all the files inside Data folder
info = struct2cell(dir('ECG_Seteb_extracted/*.mat'));

% Creat a list with the names of the files
% Currently 60 files, 6 conditions (PRE/POST x EROTIC/NEUTRAL/RESTING) for 10 participants
names = info(1,:);

NSubjects = size(names);
NSubjects = NSubjects(2)/6; % Divide by 6 conditions to get total number of participants


%% Sampling frequency 

Fs = 1000;

%% Filter parameters

% Band pass filter ---------------------------------

% Order 
n = 50;

% Low pass for the bandpass filter
lp = 1.5;

% High pass for the bandpass filter
hp = 95;

% Coefficients of filter
b = fir1(n, [lp/Fs , hp/Fs],'bandpass'); %bandpass filter with the already set values n, lp and hp


% Moving Average Filter -----------------------------

% Order
n1 = 10;

% Coefficients
b1 = (1/n1)*ones(n1,1); 


%% Load file

% Plotting flag (if 1 plots, othewise no plots)
p_flag = 0;

% Saved (arrays to save the information)
hr_all = [];
hrv_all = [];
t_all = [];
group = [];

for i = 1:length(names)

    
        clc;
        fprintf('File %i/%i \n',i,length(names)) % Counter for the files 

        % Pre-processing ------------------------------------------------------

       % Temporal file name
        file = names{i};
    
        % Read data with more double precision (double())
        ECG = double(struct2array(load(file)));
    
        % Time vector
        t = (1:1:length(ECG))*(1/Fs);
        
    
        % Get frequency characteristics
    
        % Power spectrum
        F = abs(fft(ECG));
        F = F(1:round(end/2));
    
        % Frequnecy domain 
        f = (1:1:length(F))*((Fs/2)/length(F)); 
    
        % Filter data in 100 Hz
        ECGf = filtfilt(b,1, ECG); %filtfilt goes on both direction, giving 0 phase distortion

        % Remove LF components (using a moving average filter)
        ECGf = ECGf - filtfilt(b1,1,ECGf);
    
        %Z score 
        ECGf = zscore(ECGf); 
        
            
        % Processing ----------------------------------------------------------
    
        % Findpeaks for getting HRV
    
        % Parameters for finding peaks
    
        % Minimum peak height
        mph = 4*mean(abs(ECGf)); 
        %any values can be used, depends on the data 
    
        % Minimum peak distance (based on physiological behaviour)
        mpd = (1/3)*Fs; 
    
        % Find R peaks
        [pks,locs] = findpeaks(ECGf,'MinPeakHeight',mph,'MinPeakDistance',mpd);   
    
        %calculate HR
        beats = size(pks); %dimentions
        beats = beats(1); %number of peaks
        t_signal = ((length(ECGf))/Fs)/60;
        hr = beats/t_signal; %beats per minute
        hr_all = [hr_all hr];
    
        % Get HRV  
        hrv = diff(t(locs));
    
        % Interpoalte HRV
    
        % Frequnecy of interpolation
        Fs_int = 10; %Hz  
    
        % Remove outliers
    
        % Locations of the outliers
        [~, p] =  rmoutliers(hrv,'median');
        hrv_t = (1./hrv)*60;
        
        % Remove
        hrv = hrv(~p);
    
        % Remove from time the outliers
        t1 = t(locs(1:end-1));
        t1 = t1(~p);
    
        % Time vector for interpolation
        t_hrv = 1:1/Fs_int:t(end);
        
        % Interpolation of HRV in t_hrv points       
        hrv_int = interp1(t1,hrv,t_hrv,'spline'); % Interpolate the peaks with spline method
        hrv_int = (1./hrv_int)*60; % calculate HR
    
        % Location of peaks in hrv_int  
        locs_hrv = interp1(t_hrv,hrv_int,t(locs)); % Interpolate the location

        % Crop signal for correlation, new total length is 90s instead of 120s
        lim1 = 10;
        lim2 = 100;
        
        % Cropped signal
        t_hrv_crop = t_hrv(t_hrv > lim1 & t_hrv < lim2);
        hrv_int_crop = hrv_int(t_hrv > lim1 & t_hrv < lim2);
        
        % Save croped signal
        hrv_all = [hrv_all ; hrv_int_crop];
        t_all = [t_all ; t_hrv_crop];
        
        % Save group of the signal (1:RESTING / 2:EROTIC / 3:NEUTRAL)
        
        if contains(file,'resting')           
            group = [group 1];            
        elseif  contains(file,'erotico')            
            group = [group 2];
        elseif contains(file,'neutro')
            group = [group 3];
        end
        
        
        % Plot         
        if p_flag == 1 % If 1 it means plots are plotted 

            figure

            subplot(3,2,[1 2])

            hold on    
            plot(t,ECGf)
            xlabel('Time (s)')
            ylabel('Amplitude')
            scatter(t(locs),ECGf(locs),20,'filled')
            yline(mph,'--r')
            title('Filtered Data') 
            legend('ECG','R Peak Detected','Min Peak Height','Location','southeast')
            xlim([0 20])
            hold off

            subplot(3,2,[3 4])
            hold on
            plot(t_hrv, hrv_int)
            plot(t(locs(1:end-1)),hrv_t)
            scatter(t(locs),locs_hrv,20,'filled')
            legend('Outliers removed','Normal','R Peaks')
            hold off
            xlim([0 20])
            ylim([ mean(hrv_int(~isnan(hrv_int)))-10, mean(hrv_int(~isnan(hrv_int)))+10])
            title('HRV')

            subplot(3,2,5)
            plot(t,ECG)
            xlabel('Time (s)')
            ylabel('Amplitude')
            title('Raw Data')
            xlim([0 20])

            subplot(3,2,6)
            plot(f,log(F))
            xlabel('Frequency Hz')
            ylabel('Magnitude')
            xlim([0 110])
            title('Power Raw Signal')

            sgtitle(file) %(i) for file place in the list
        end
       

end

hr_all = hr_all';

% Save the hr and hrv matrix and group vector
filename_hr = 'D:\Users\PERSONAL\Documents\psychologie Master thesis\ECG_Seteb_extracted\HR_all_Pt.mat';
save(filename_hr, 'hr_all');

filename_hrv = 'D:\Users\PERSONAL\Documents\psychologie Master thesis\ECG_Seteb_extracted\HRV_all_Pt.mat';
save(filename_hrv,'hrv_all')

filename_hr = 'D:\Users\PERSONAL\Documents\psychologie Master thesis\ECG_Seteb_extracted\Group_Pt.mat';
save(filename_hr, 'group');

%% Correlation 

ISC_all = [];

% Select by group and get correlation  matrix
groups = ["Resting","Erotico","Neutro"];

for j = 1:2 %for the sessions

    figure;
    plot_count = 0;

    for i = 1:3 %for the conditions

        % Get matrix of HRV per group
        tmp_hrv = hrv_all(group == i,:)'; %Get all participants of 1 condition
        tmp_hrv = tmp_hrv(:,j:2:end); % select only the files of one of the 
        % time points. start either at the first or second file and then 
        % select one file, one not, one file, one not. This is the order in 
        % which the files are in the folder. 
        % eg. erotic_T1, erotic_T2, neutral_T1, neutral_

        % Get correlation
        tmp_corr_mat = atanh(corrcoef(tmp_hrv)); 
        s = size(tmp_corr_mat);
        index_diag = 1:s(1)+1:s(1)*s(2);
        tmp_corr_mat(index_diag) = NaN; % Set diagonal to Nan

        % Calculate the average correlation per participant, autocorrelation is not included
        ISC = tanh(nanmean(tmp_corr_mat));
        ISC_all = [ISC_all ISC]; % Stores first moment 1 (resting -> erotic -> neutral) and then the same for moment 2

        if i > 1 %Just plot erotic and neutral
            plot_count = plot_count + 1;    
            subplot(2,2,plot_count)
           
            % Show matrix
    
            imagesc(tmp_corr_mat)
            colorbar
            ylabel('Subject no.')
            xlabel('Subject no.')
            title(groups(i))
    
            plot_count = plot_count + 1;
    
            subplot(2,2,plot_count)
    
            % Plot
            stem(ISC,'filled') 
            xlabel('Subject no.')
            ylabel('ISC-HR')
            title(groups(i))
        end 
    end

end

% Save the ISC matrix 
filename_ISC = 'D:\Users\PERSONAL\Documents\psychologie Master thesis\ECG_Seteb_extracted\ISC_all_Pt.mat';
save(filename_ISC,'ISC_all');

%% Graph HR
% Show only descriptive statistics. Compare patients to controls for the
% different conditions. 
% Maybe T-test between neutral an erotic for healthy sample?

% Separate HR by condition
%resting_hr  = hr_all(group == 1); 
erotic_hr   = hr_all(group == 2);
neutral_hr  = hr_all(group == 3);

% Separate HR by session 
erotic_hr_ses1  = erotic_hr(1:2:end)
erotic_hr_ses2  = erotic_hr(2:2:end)
neutral_hr_ses1 = neutral_hr(1:2:end)
neutral_hr_ses2 = neutral_hr(2:2:end)

%% Graph HRV
% Separate HRV by condition
%resting_hrv = hrv_all(group == 1,:)'; 
erotic_hrv  = hrv_all(group == 2,:)';
neutral_hrv = hrv_all(group == 3,:)';

% Calculate the sd of hr based on the hrv
erotic_sd_hr  = std(erotic_hrv,0,2); %by row, if 1 then by column
neutral_sd_hr = std(neutral_hrv,0,2); %by row

% Separate by session 
erotic_sd_hr_ses1  = erotic_sd_hr(1:2:end) 
erotic_sd_hr_ses2  = erotic_sd_hr(2:2:end) 
neutral_sd_hr_ses1 = neutral_sd_hr(1:2:end) 
neutral_sd_hr_ses2 = neutral_sd_hr(2:2:end) 

%% Statistics ISC 
% Separate ISC by condition
% resting_ISC  = ISC_all(group == 1);
erotic_ISC   = ISC_all(group == 2);
neutral_ISC  = ISC_all(group == 3);

% Separate ISC by session 
erotic_ISC_ses1  = erotic_ISC(1:2:end)
erotic_ISC_ses2  = erotic_ISC(2:2:end)
neutral_ISC_ses1 = neutral_ISC(1:2:end)
neutral_ISC_ses2 = neutral_ISC(2:2:end)

% Significance Test ISC
Npermutations = 10; %10000
hrv_shifted = hrv_all;
T = size(hrv_all,2);


% Compute ISC based on random circular shifts
for iRandShift=Npermutations:-1:1 
          
    % Correlation 
    % Select by group and get correlation  matrix
    
    ISC_partial_shifted = [];
       
    for j = 1:2 %for sessions 

        for i = 2:3% for de    
    
            % Get matrix of HRV per group
            tmp_hrv = hrv_shifted(group == i,:)'; 
            tmp_hrv = tmp_hrv(:,j:2:end);

            %circular shuffle, no tricks. 
            for iSubject=1:NSubjects
                tmp_hrv(iSubject,:) = circshift(tmp_hrv(iSubject,:),round(T*rand(1))); % Give every participant a different shift
            end
             
            % Get correlation
            tmp_corr_mat = atanh(corrcoef(tmp_hrv)); 
            s = size(tmp_corr_mat);
            index_diag = 1:s(1)+1:s(1)*s(2);
            tmp_corr_mat(index_diag) = NaN; % Set diagonal to Nan

            % Calculate the average correlation per participant, autocorrelation is not included
            ISC = tanh(nanmean(tmp_corr_mat));
            ISC_partial_shifted = [ISC_partial_shifted ISC]; % Stores first moment 1 (resting -> erotic -> neutral) and then the same for moment 2
        
        end
   
    end
    
    ISC_all_shifted(:,iRandShift)= ISC_partial_shifted;  

end

% Separate ISC shifted by condition
% resting_ISC_shift  = ISC_all(group == 1);
erotic_ISC_shift   = ISC_all_shifted(group == 2);
neutral_ISC_shift  = ISC_all_shifted(group == 3);

% Separate ISC by session 
erotic_ISC_shift_ses1  = erotic_ISC_shift(1:2:end)
erotic_ISC_shift_ses2  = erotic_ISC_shift(2:2:end)
neutral_ISC_shift_ses1 = neutral_ISC_shift(1:2:end)
neutral_ISC_shift_ses2 = neutral_ISC_shift(2:2:end)

% Percentiles based on the random circular shift distribution 
    % Erotic
erotic_P95_ses1 = prctile(erotic_ISC_shift_ses1,95,2);
erotic_P5_ses1  = prctile(erotic_ISC_shift_ses1,5,2);
erotic_P95_ses2 = prctile(erotic_ISC_shift_ses2,95,2);
erotic_P5_ses2  = prctile(erotic_ISC_shift_ses2,5,2);
    % Neutral 
neutral_P95_ses1 = prctile(neutral_ISC_shift_ses1,95,2);
neutral_P5_ses1  = prctile(neutral_ISC_shift_ses1,5,2); 
neutral_P95_ses2 = prctile(neutral_ISC_shift_ses2,95,2);
neutral_P5_ses2  = prctile(neutral_ISC_shift_ses2,5,2);

% Calculate the p-value 
pval_pos_neutral_ses1 = [];
pval_neg_neutral_ses1 = [];
pval_pos_erotic_ses1 = [];
pval_neg_erotic_ses1 = [];
pval_pos_neutral_ses2 = [];
pval_neg_neutral_ses2 = [];
pval_pos_erotic_ses2 = [];
pval_neg_erotic_ses2 = [];


% P value
    % neutral session 1
for iData=1:length(neutral_ISC_shift_ses1)
   
    pval_pos_neutral_ses1(iData) = sum(mean(neutral_ISC_ses1(iData)<neutral_ISC_shift_ses1(iData,1:Npermutations))); %compare (<) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
    pval_neg_neutral_ses1(iData) = sum(mean(neutral_ISC_ses1(iData)>neutral_ISC_shift_ses1(iData,1:Npermutations))); %compare (>) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
       
end 
pval_pos_neutral_ses1 = pval_pos_neutral_ses1/Npermutations; %percentiel 95
pval_neg_neutral_ses1 = pval_neg_neutral_ses1/Npermutations; %percentiel 5


    % erotic session 1
for iData=1:length(erotic_ISC_shift_ses1)
   
    pval_pos_erotic_ses1(iData) = sum(mean(erotic_ISC_ses1(iData)<erotic_ISC_shift_ses1(iData,1:Npermutations))); %compare (<) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
    pval_neg_erotic_ses1(iData) = sum(mean(erotic_ISC_ses1(iData)>erotic_ISC_shift_ses1(iData,1:Npermutations))); %compare (>) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
       
end
pval_pos_erotic_ses1 = pval_pos_erotic_ses1/Npermutations; %percentiel 95
pval_neg_erotic_ses1 = pval_neg_erotic_ses1/Npermutations; %percentiel 5 


    % neutral session 2
for iData=1:length(neutral_ISC_shift_ses2)
   
    pval_pos_neutral_ses2(iData) = sum(mean(neutral_ISC_ses1(iData)<neutral_ISC_shift_ses2(iData,1:Npermutations))); %compare (<) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
    pval_neg_neutral_ses2(iData) = sum(mean(neutral_ISC_ses1(iData)>neutral_ISC_shift_ses2(iData,1:Npermutations))); %compare (>) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
       
end 
pval_pos_neutral_ses2 = pval_pos_neutral_ses2/Npermutations; %percentiel 95
pval_neg_neutral_ses2 = pval_neg_neutral_ses2/Npermutations; %percentiel 5

   
    % erotic session 2
for iData=1:length(erotic_ISC_shift_ses2)
   
    pval_pos_erotic_ses2(iData) = sum(mean(erotic_ISC_ses2(iData)<erotic_ISC_shift_ses2(iData,1:Npermutations))); %compare (<) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
    pval_neg_erotic_ses2(iData) = sum(mean(erotic_ISC_ses2(iData)>erotic_ISC_shift_ses2(iData,1:Npermutations))); %compare (>) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
       
end
% P value
pval_pos_erotic_ses2 = pval_pos_erotic_ses2/Npermutations; %percentiel 95
pval_neg_erotic_ses2 = pval_neg_erotic_ses2/Npermutations; %percentiel 5 

%% Graph ISC 









