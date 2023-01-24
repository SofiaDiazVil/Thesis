clc;close all;clear

%% Add to path Data folder (data experiment and time stamps)
addpath('./experiment_gent/ECG_exp/');
addpath('./experiment_gent/TimeStamps_exp/files/');


%% Read files
% Read name of all the files inside folder
    %ECG data
ECG_files  = struct2cell(dir('experiment_gent/ECG_exp/*.txt')); 
    %timestamps
timestamps = dir('experiment_gent/TimeStamps_exp/files/*'); 
timestamps=natsortfiles(timestamps); %sort the files in the correct order 
timestamps = struct2cell(timestamps(3:32,:)); %exclude other files 

% Creat a list with the names of the files
names_ecg  = ECG_files(1,:);
names_time = timestamps(1,:);

%% Relevant Variables 
% Get the number of subjects 
NSubjects = size(names_ecg); %30 
% Set the sample frequency
Fs = 10000; 


%% Cut the data into the correct segments and save the files 
% Create time vector (10 minute signal)
ECG_sample = double(readmatrix(names_ecg{1}));
t = (1:1:length(ECG_sample))*(1/Fs);

for i = 1:length(names_time)

    Nsub = i+1; %participant 1 was lost, so data starts from participant 2
    name_file_resting = strcat('resting_ecg_subj', string(Nsub),'.mat');
    name_file_erotic  = strcat('erotic_ecg_subj', string(Nsub),'.mat');
    name_file_neutral = strcat('neutral_ecg_subj', string(Nsub),'.mat');

   %import the files
    clc; 
    file_time = names_time{i}; 
    file_ECG  = names_ecg{i};
    timestamp = readmatrix(file_time);
    ECG = double(readmatrix(file_ECG));

    % Resting state 
    if timestamp(1,2) == -99 %everyone has resting as the first condition
        resting_times =  timestamp(1,3:4);
   
    end 
    
    % Erotic condition 
    if timestamp(2,2) == 0 
        erotic_times = timestamp(2,3:4);
    
    elseif timestamp(3,2) == 0
        erotic_times = timestamp(3,3:4);
    
    end
        
    % Neutral condition 
    if timestamp(2,2) == 1
        neutral_times = timestamp(2,3:4);
    
    elseif timestamp(3,2) == 1
        neutral_times = timestamp(3,3:4);
    
    end 
 
    % Cropped signals
        
        % resting
    cropped_ECG_resting = ECG(t > resting_times(1) & t < resting_times(2));
    save(name_file_resting,'cropped_ECG_resting');
        
        % erotic 
    cropped_ECG_erotic = ECG(t > erotic_times(1) & t < erotic_times(2));
    save(name_file_erotic,'cropped_ECG_erotic');
        
        % neutral 
    cropped_ECG_neutral = ECG(t > neutral_times(1) & t < neutral_times(2));
    save(name_file_neutral,'cropped_ECG_neutral');
      
end 

%% Preprocessing preparation 

% Load cropped files 
addpath('D:\Users\PERSONAL\Documents\psychologie Master thesis\experiment_gent\Cropped ECG Signals');
Data  = struct2cell(dir('D:\Users\PERSONAL\Documents\psychologie Master thesis\experiment_gent\Cropped ECG Signals\*.mat')); 

% Create a list with the names of the files
names = Data(1,:);

% Set the sample frequency (when running from here)
Fs = 10000; 

% Filter parameters 

% Band pass filter --------------------
% Order 
n = 50;
% Low pass for the bandpass filter
lp = 1.5;
% High pass for the bandpass filter
hp = 95;
% Coefficients of filter
b = fir1(n, [lp/Fs , hp/Fs],'bandpass'); %bandpass filter with the already set values n, lp and hp

% Notch filter 50 Hz -------------------- 
f_line=50;
fcomb = [[45 49 51 54], [45 49 51 54]+f_line, [45 49 51 54]+2*f_line [45 49 51 54]+3*f_line [45 49 51 54]+4*f_line];
mags = [[1 0 1], [0 1], [0 1],  [0 1], [0 1]];
dev = [[0.5 0.1 0.5], [0.1 0.5], [0.1 0.5], [0.1 0.5], [0.1 0.5]];
[n,Wn,beta,ftype] = kaiserord(fcomb,mags,dev,Fs);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');


%% Preprocessing and Processing 

% Plotting flag (if 1 plots, othewise no plots)
p_flag = 0;

% Saved (arrays to save the information)
hr_all = [];
hrv_all = [];
t_all = [];
group = [];

% Loop for the processing 
for i = 1:length(names) %repeats for all files in the folder 
    
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

    %Notch filter in 50 Hz
    ECGf=filtfilt(hh,1,ECGf);

    %Z score 
    ECGf = zscore(ECGf); 
    
    % Additional filtering is needed 

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
    hrv = diff(t(locs)); %time between peaks

    % Interpolate HRV

    % Frequnecy of interpolation
    Fs_int = 10; %Hz  

    % Remove outliers

    % Locations of the outliers
    [~, p] =  rmoutliers(hrv); % using default median method, returns data without outliers 
    hrv_t = (1./hrv)*60; % calculate HR
    
    % Remove
    hrv = hrv(~p); % exclude outliers from data 

    % Remove from time the outliers
    t1 = t(locs(1:end-1));
    t1 = t1(~p);

    % Time vector for interpolation
    t_hrv = 1:1/Fs_int:t(end);
    
    % Interpolation of HRV in t_hrv points       
    hrv_int = interp1(t1,hrv,t_hrv,'spline'); % Interpolate the peaks with spline method
    hrv_int = (1./hrv_int)*60; % calculate HR

    % Location of peaks in hrv_int  
    locs_hrv = interp1(t_hrv,hrv_int,t(locs)); % Interpolate the location/time

    % Crop signal for correlation, new total length is 95s 
    lim1 = 5;
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
    elseif  contains(file,'erotic')            
        group = [group 2];
    elseif contains(file,'neutral')
        group = [group 3];
    end

    % Plots         
    if p_flag == 1 % If 1 plots are plotted 

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

% Save the hr and hrv matrix 
filename_hr = 'D:\Users\PERSONAL\Documents\psychologie Master thesis\experiment_gent\HR_ALL.mat';
save(filename_hr, 'hr_all');

filename_hrv = 'D:\Users\PERSONAL\Documents\psychologie Master thesis\experiment_gent\HRV_ALL.mat';
save(filename_hrv,'hrv_all')


%% Correlation 

ISC_all = [];

% Select by group and get correlation  matrix
groups = ["resting","erotic","neutral"];

figure;
plot_count = 0;

for i = 1:3    

    % Get matrix of HRV per group
    tmp_hrv = hrv_all(group == i,:)'; %Get all participants of 1 condition

    % Get correlation
    tmp_corr_mat = atanh(corrcoef(tmp_hrv)); 
    s = size(tmp_corr_mat);
    index_diag = 1:s(1)+1:s(1)*s(2);
    tmp_corr_mat(index_diag) = NaN; % Set diagonal to Nan
    
    plot_count = plot_count + 1;    
    subplot(3,2,plot_count)
   
    % Show matrix

    imagesc(tmp_corr_mat)
    colorbar
    ylabel('Subject no.')
    xlabel('Subject no.')
    title(groups(i))

    plot_count = plot_count + 1;

    subplot(3,2,plot_count)

    % Plot
     % Calculate the average correlation per participant, autocorrelation is not included
    ISC = tanh(nanmean(tmp_corr_mat));
    ISC_all = [ISC_all ISC]; % Stores each condition (resting -> erotic -> neutral)
    stem(ISC,'filled') 
    xlabel('Subject no.')
    ylabel('ISC-HR')
    title(groups(i))
end

% Save the ISC matrix 
filename_ISC = 'D:\Users\PERSONAL\Documents\psychologie Master thesis\experiment_gent\ISC_ALL.mat';
save(filename_ISC,'ISC_all');


%% Statistics HR
% Show only descriptive statistics. Compare patients to controls for the
% different conditions. 
% Maybe T-test between neutral an erotic for healthy sample?

% Separate HR by condition
resting_hr  = hr_all(group == 1);
erotic_hr   = hr_all(group == 2);
neutral_hr  = hr_all(group == 3);

% Two tailed t-test 
[h_hr,p_hr,ci_hr,stats_hr] = ttest(erotic_hr,neutral_hr,'Alpha',0.1); % erotic-neutral


%% Statistics HRV
% Show only descriptive statistics. Look at the standard deviation of HR 
% an compare between patients andcontrols for the different conditions 
% Maybe T-test between neutral an erotic for healthy sample?

% Separate HRV by condition
resting_hrv = hrv_all(group == 1,:)';
erotic_hrv  = hrv_all(group == 2,:)';
neutral_hrv = hrv_all(group == 3,:)';

%calculate the sd of hr based on the hrv
erotic_sd_hr  = std(erotic_hrv,0,2); %by row, if 1 then by column
neutral_sd_hr = std(neutral_hrv,0,2); %by row

% Two tailed t-test 
[h_hrv,p_hrv,ci_hrv,stats_hrv] = ttest(erotic_sd_hr,neutral_sd_hr,'Alpha',0.1); % erotic-neutral


%% Statistics ISC 
% Separate ISC by condition
resting_ISC  = ISC_all(group == 1);
erotic_ISC   = ISC_all(group == 2);
neutral_ISC  = ISC_all(group == 3);

% Make a paired distribution plot 
    % Select the data
%Sample  = 1+zeros(1,30); % 1 is healthy 0 is patient 
Dataset = [erotic_ISC' neutral_ISC']; %dataset to graph

    % Function
%coorddata = [1 2]; %matrix with the variables you want to graph
p = parallelplot(Dataset);
% CoordinateData = select the bars/groups and the order they will be presented in
% GroupData = variable by wich you will group the data (color lines)
p.CoordinateTickLabels = {'Erotic','Neutral'}; %assign the label of the bars


%% Significance Test ISC
Npermutations = 10; %10000
NSubjects     = 30
hrv_shifted = hrv_all;
T = size(hrv_all,2);


% Compute ISC based on random circular shifts
for iRandShift=Npermutations:-1:1 
          
    % Correlation 
    % Select by group and get correlation  matrix
    
    ISC_partial_shifted = [];
       
      for i = 1:3    

        % Get matrix of HRV per group
        tmp_hrv = hrv_shifted(group == i,:)'; 
        
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
    
    ISC_all_shifted(:,iRandShift)= ISC_partial_shifted;  % Store the shifted ISC for each repetition

end

% Calculate the percentiles of the distribution 
% interval: calculate the percentiles (5 and 95 would be fine) of the
% distribution and compare your values to that 
P95 = prctile(ISC_all_shifted,95,2);
P5  = prctile(ISC_all_shifted,5,2);
%%
% Graph the distribution 
% Select the data
% Dataset = [erotic_ISC neutral_ISC]; %data set to graph
% Cond1   = 1+zeros(1,30); %erotic
% Cond2   = 2+zeros(1,30); %neutral
% CondAll = [Cond1 Cond2];
% Dataset = [CondAll; Dataset];
% Violin graph
Violin_graph = Violin(Data, Position); % Data might be the dataset or maybe the groups? position is a vector with position
%Violin_graph.ScatterPlot2 = ;%plot 2 sets of data 
Violin_graph.ShowMean = 'true';
Violin_graph.HalfViolin = 'left'
Violin_graph.ShowData = 'false'
%%

% Calculate the p-value 
pval_pos = [];
pval_neg = [];
% Loop to compare 
for iData=1:length(ISC_all_shifted)
   
    pval_pos(iData) = sum(mean(ISC_all(iData)<ISC_all_shifted(iData,1:Npermutations))); %compare (<) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
    pval_neg(iData) = sum(mean(ISC_all(iData)>ISC_all_shifted(iData,1:Npermutations))); %compare (>) ISC with shifted ISC (1 or 0), add them up to know how many fit the condition
       
end 

% P value: for a non parametric distribution you count the amount of values
% higher than your value in the distribution and divide that by the total 
% amount of values in the distribution. 
pval_pos = pval_pos/Npermutations; %percentiel 95
pval_neg = pval_neg/Npermutations; %percentiel 5









