%% Insula Data Processing Script %%
% This code loads and analyzes data from the Resting State recordings.
%
% Written by John Myers

%% Open file and Establish Dimensionality of Input Data
% load('CAR_insula.mat'); 
% hippocampal_channels = [1:20];
insula_channels = [143	144	174	175	176	177	178	179	180	181	182	183	184]; % insula channels
srate = 2000; % ecog sampling rate
ecog_file_ID = 'YCQDatafile2031.ns3'; %
ECoG_data = openNSx(ecog_file_ID, 'uV'); % read in EMU the data
data = double(ECoG_data.Data); % 1 NSP

% NSPs
% data = double(ECoG_data.Data{1,2}); 
% data = [data1(:, 1:length(data1)); data(:, 1:length(data1))]; % combine
avg_data = mean(data);
insula_data = data(insula_channels, :); % select insula channels

% temporal parameters
data_length = size(insula_data, 2); 
data_length_secs = data_length/srate; 
data_length_minutes = data_length_secs/60; 

%% Set Up Probes and Re-reference (BIPOLAR REFERENCE)
probe_channels = {1:11}; % create cell array where each row defines the channels for each probe
% probe_channels = {1:9; 10:20}; % create cell array where each row defines the channels for each probe
reref_data     = bipolar_ref(data, probe_channels); % re-reference all the probes and create structure

%% Set Up Probes and Re-reference (COMMON AVERAGE REFERENCE)
CAR_data = insula_data - avg_data;

%% Filter Out Line Noise and Check Spectrum
probe_number = 1; 
channel = 1; 
probe_data = reref_data(probe_number).probe(:, :); 
% probe_data = CAR_data;

notch_center = 60;
disp(['Notch Filtering Out ', num2str(notch_center), ' Hz...']); 
wo = notch_center/(srate/2); % notch center (radians/second)
bw = wo/30; % bandwidth
[b, a] = iirnotch(wo, bw); % notch filter design
% notched_data = filtfilt(b, a, double(probe_data)');
notched_data = filtfilt(b, a, probe_data');

% notch out 120 Hz
notch_center2 = 120;
disp(['Notch Filtering Out ', num2str(notch_center2), ' Hz...']); 
wo = notch_center2/(srate/2); % notch center (radians/second)
bw = wo/50; % bandwidth
[b, a] = iirnotch(wo, bw); % notch filter design
notched_data = filtfilt(b, a, notched_data);

% notch out
notch_center = 56;
disp(['Notch Filtering Out ', num2str(notch_center), ' Hz...']); 
wo = notch_center/(srate/2); % notch center (radians/second)
bw = wo/50; % bandwidth
[b, a] = iirnotch(wo, bw); % notch filter design
notched_data = filtfilt(b, a, notched_data);

% notch out
notch_center = 63;
disp(['Notch Filtering Out ', num2str(notch_center), ' Hz...']); 
wo = notch_center/(srate/2); % notch center (radians/second)
bw = wo/50; % bandwidth
[b, a] = iirnotch(wo, bw); % notch filter design
% notched_data = filtfilt(b, a, double(probe_data)');
notched_data = filtfilt(b, a, notched_data);

probe_data = notched_data';
figure; pspectrum(probe_data', srate, 'FrequencyLimits', [1, 150]);
% reref_data(probe_number).probe(:, :) = probe_data; % save data in structure

%% Concatenate Data (BIPOLAR)
nprobes = 1; 
for p = 1:nprobes
    dat = reref_data(p).probe; 
    if p == 1
       referenced_data = dat;
    else        
       referenced_data = [referenced_data; dat];
    end
end

% Save Referenced Data
YDH_data = referenced_data;

%% Store CAR data (CAR)
coordinates = [];
subj = 9;
CAR_insula(subj).data = probe_data;

%%
CAR_insula(subj).coordinates = coordinates;

%% Append Datasets ACROSS SUBJECTS and Find Mean Frequency Bands (MODAL)
% % Left Hemisphere
load YCS_data2;     load YCS_eloc_MNI;
load YCQ_data_left; load YCQ_eloc_MNI_left;
load YCY_data_left; load YCY_eloc_MNI_left;
load YDD_data;      load YDD_eloc_MNI; 
load YDF_data_left; load YDF_eloc_MNI_left; 
load YDJ_data;      load YDJ_eloc_MNI; 
load YDH_data_left; load YDH_eloc_MNI_left; 

insula(1).data = YCS_data; 
insula(2).data = YCQ_data_left; 
insula(3).data = YCY_data_left; 
insula(4).data = YDD_data; 
insula(5).data = YDF_data_left; 
insula(6).data = YDH_data_left; 
insula(7).data = YDJ_data; 

% Insula Contacts and Time Lengths of Each Recording
contacts = [size(YCS_data, 1); size(YCQ_data_left, 1); size(YCY_data_left, 1); size(YDD_data, 1); size(YDH_data_left, 1); size(YDF_data_left, 1); size(YDJ_data, 1)];
time_lengths = [size(YCS_data, 2); size(YCQ_data_left, 2); size(YCY_data_left, 2); size(YDD_data, 2); size(YDF_data_left, 2); size(YDH_data_left, 2); size(YDJ_data, 2)]; 

% Coordinates
eloc = [YCS_eloc; YCQ_eloc_left; YCY_eloc_left; YDD_eloc; YDF_eloc_left; YDH_eloc_left; YDJ_eloc];

% Right Hemisphere
% load YCP_data; load YCP_eloc_MNI; 
% load YCQ_data; load YCQ_eloc_MNI; 
% load YCU_data; load YCU_eloc_MNI; 
% load YCW_data; load YCW_eloc_MNI; 
% load YCY_data; load YCY_eloc_MNI; 
% load YDF_data; load YDF_eloc_MNI; 
% load YDH_data; load YDH_eloc_MNI; 
% 
% insula(1).data = YCP_data; 
% insula(2).data = YCQ_data; 
% insula(3).data = YCU_data; 
% insula(4).data = YCW_data; 
% insula(5).data = YCY_data; 
% insula(6).data = YDF_data; 
% insula(7).data = YDH_data; 

% Insula Contacts and Time Lengths of Each Recording
% contacts = [size(YCP_data, 1); size(YCQ_data, 1); size(YCU_data, 1); size(YCW_data, 1); size(YCY_data, 1); size(YDF_data, 1); size(YDH_data, 1)]; 
% time_lengths = [size(YCP_data, 2); size(YCQ_data, 2); size(YCU_data, 2); size(YCW_data, 2); size(YCY_data, 2); size(YDF_data, 2); size(YDH_data, 2);];

% Coordinates
% eloc = [YCP_eloc; YCQ_eloc; YCU_eloc; YCW_eloc; YCY_eloc; YDF_eloc; YDH_eloc;];

% Append Data
n_subj   = size(contacts, 1); 
min_time = min(time_lengths);
n_contacts = sum(contacts);

for c = 1:n_subj
    contact_data = insula(c).data(:, 1:min_time);
    if  c == 1
        data = contact_data;
    else 
        data = [data; contact_data];
    end
end

% Time-freq Parameters
srate = 2000;
freq_range = [1 50];
params.srate = srate; % sampling rate
params.wavefreqs = freq_range(1):0.1:freq_range(2); % frequency
params.wavecycles = 6; % wavelet cycles
params.local_winsize_sec = 10; % window size (secs)

% MODAL for Freq Range
clear bands
[frequency_sliding, bands, bandpow, bandphases, pow] = MODAL(mean(data), params); 

%% MODAL - LOAD Data and Find Freq Bands (V1)
% load('YDF_data_left.mat'); 
% load('YDF_eloc_MNI_left.mat'); 
% eloc = YDF_eloc_left;

% load('INSULA_ALL_LH.mat', 'INSULA_ALL_LH'); 
% load('INSULA_ALL_LH.mat', 'INSULA_ALL_LH_coords');
% load('INSULA_ALL_RH.mat', 'INSULA_ALL_RH'); 
% load('INSULA_ALL_RH.mat', 'INSULA_ALL_RH_coords');
% 
% raw_data = INSULA_ALL_RH;
% eloc = INSULA_ALL_RH_coords;
% x_coord  = eloc(:,1); 
% y_coord  = eloc(:,2); 
% z_coord  = eloc(:,3); 
% 
% srate = 2000;
% params.srate = srate; % sampling rate
% params.wavefreqs = 1:0.1:50; % frequency
% params.wavecycles = 6; % wavelet cycles
% params.local_winsize_sec = 10; % window size (secs)
% 
% % MODAL for Freq Range
% clear bands
% [frequency_sliding, bands, bandpow, bandphases, pow] = MODAL(mean(raw_data), params); 

%% MODAL - LOAD Data and Find Freq Bands (based on median across both hemispheres)
load('INSULA_ALL_LH.mat', 'INSULA_ALL_LH'); 
load('INSULA_ALL_RH.mat', 'INSULA_ALL_RH'); 
raw_data = [INSULA_ALL_LH; INSULA_ALL_RH];

% parameters
srate = 2000;
params.srate = srate; % sampling rate
params.wavefreqs = 1:0.1:50; % frequency
params.wavecycles = 6; % wavelet cycles
params.local_winsize_sec = 10; % window size (secs)

% Run MODAL
[~, bands, ~, ~, pow] = MODAL(median(raw_data), params); 
logPow   = log(mean(pow, 2)); 
logFreqs = log(params.wavefreqs);
nbands = size(bands,1);

% plot
figure;
psd_sem  = std(logPow, 1)/sqrt(size(raw_data, 1));
[b, output] = robustfit(logFreqs, logPow');
fit_line = b(1) + b(2) .* logFreqs;
h = boundedline(double(params.wavefreqs'), double(logPow), double(psd_sem), 'cmap', 'b', 'alpha');
h.LineWidth = 3.0; 
hold on; 
plot(params.wavefreqs, fit_line); 

bands = [6,9; 15,32];

for band = 1:nbands
    xpoints = [bands(band, 1), bands(band, 1), bands(band, 2), bands(band, 2)];
    ypoints = [0, max(logPow), max(logPow), 0];
    hold on; 
    a = fill(xpoints, ypoints, 'k');
    a.FaceAlpha = 0.05;
end        
set(gcf, 'color', [1 1 1]);
set(gca, 'color', [1 1 1]);
set(gca, 'fontsize', 35);
xlim([1 50]);

%% SPECTRAL PEAK ANALYSIS
clear all_power input_data peak_frequencies peak_power;
raw_data = INSULA_ALL_RH; % insert variable as input data
eloc     = INSULA_ALL_RH_coords; % inser variable as electrode coordinates (RAS coordinate system)
x_coord  = eloc(:,1); 
y_coord  = eloc(:,2); 
z_coord  = eloc(:,3); 

freq_range = bands(2,:); 
params.wavefreqs = freq_range(1):0.1:freq_range(2); % frequency
nchannels = size(raw_data, 1); 
color_set = distinguishable_colors(nchannels,'k'); 

% compute PSD
figure; 
for chan = 1:nchannels
    disp(['Computing PSD from ', num2str(freq_range(1)), '-', num2str(freq_range(2)), ' Hz for channel ', num2str(chan), ' ...']);
    [phase, pow] = multiphasevec2(params.wavefreqs, raw_data(chan, :), params.srate, params.wavecycles);
    logPow   = log(mean(pow,2)); % vector of power values from wavelet averaged across time for each frequency 
    logFreqs = log(params.wavefreqs);
    [b, output] = robustfit(logFreqs, logPow');
%     fitLine = b(1) + b(2).*logFreqs; % 1/f
    [pks, locs] = findpeaks(output.resid);
    peak_freqs = params.wavefreqs(locs)'; % select peak frequencies
    [power, max_power_idx] = max(pks); % find index of maximum power
    peak_frequency         = peak_freqs(max_power_idx); % peak_frequency = frequency with greatest power    
    power_peak             = logPow(max_power_idx); % power level of peak frequency    
    
    if  ~isempty(peak_frequency)
        peak_frequencies(:, chan) = peak_frequency; % stores only peak frequencies
        peak_power(:, chan) = power_peak;
    else 
        peak_frequencies(:, chan) = NaN; % use zero if there was no peak
        peak_power(:, chan) = NaN;
    end
    all_power(chan, :) = logPow; 
    plot(params.wavefreqs, logPow, 'Color', color_set(chan,:), 'LineWidth', 5.0); hold on; 
end

% Transpose Data and Remove Missing Values
[peak_frequencies_rm, rm_idx] = rmmissing(peak_frequencies'); 
peak_power_rm = peak_power(~rm_idx)'; 
yz_coordinates = [y_coord(~rm_idx), z_coord(~rm_idx)]; 
eloc_rm = eloc(~rm_idx,:);

input_data = [peak_frequencies_rm, peak_power_rm, y_coord(~rm_idx), z_coord(~rm_idx)]; 
[corr_val, corr_pval] = corr(input_data, 'rows', 'complete');

% Run Stats on Single-Subject Data
% Create Table of Variables
var_names = {'Peak_Frequency', 'Peak_Power', 'x_coord', 'y_coord', 'z_coord'};
T = table(peak_frequencies_rm, peak_power_rm, x_coord(~rm_idx), y_coord(~rm_idx), z_coord(~rm_idx), 'VariableNames', var_names);

% Fit Models
frequency_mdl = fitlm(T, 'interactions', 'ResponseVar', 'Peak_Frequency',...
            'PredictorVars',  {'y_coord', 'z_coord'},...
            'RobustOpts', 'on'); 
        
power_mdl = fitlm(T, 'interactions', 'ResponseVar', 'Peak_Power',...
            'PredictorVars',  {'y_coord', 'z_coord'},...
            'RobustOpts', 'on');        
       
% 3D Visualizations
% FREQUENCY VISUALIZATION
figure; 
scatter3(y_coord(~rm_idx), z_coord(~rm_idx), peak_frequencies_rm, 100, 'filled'); 
% [X,Y,Z] = meshgrid(peak_frequencies, y_coord, z_coord);

% Fit a Plane and Plot
hold on; 
x = y_coord(~rm_idx); 
y = z_coord(~rm_idx); 
z = peak_frequencies_rm; 

B = [x(:) y(:) ones(size(x(:)))] \ z(:); 

xv = linspace(min(x), max(x), nchannels)';
yv = linspace(min(y), max(y), nchannels)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);
s = mesh(X, Y, Z, 'FaceAlpha', 0.5);
s.FaceColor = 'Flat';
set(gca, 'xdir', 'reverse');

xlabel('AP Axis'); 
ylabel('SI Axis'); 
zlabel('Peak Frequency (Hz)'); 
colormap(jet);
view([-20.12 34.91]);
caxis([min(z), max(z)]);
zlim(freq_range);

% POWER VISUALIZATION
figure; 
scatter3(y_coord(~rm_idx), z_coord(~rm_idx), peak_power_rm, 100, 'filled'); 
% [X,Y,Z] = meshgrid(peak_frequencies, y_coord, z_coord);

% Fit a Plane and Plot
hold on; 
x = y_coord(~rm_idx); 
y = z_coord(~rm_idx); 
z = peak_power_rm; 

B = [x(:) y(:) ones(size(x(:)))] \ z(:); 

xv = linspace(min(x), max(x), nchannels)';
yv = linspace(min(y), max(y), nchannels)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);
s = mesh(X, Y, Z, 'FaceAlpha', 0.5);
s.FaceColor = 'Flat';
set(gca, 'xdir', 'reverse');

xlabel('AP Axis'); 
ylabel('SI Axis'); 
zlabel('Peak Power (dB)'); 
colormap(jet);
view([-20.12 34.91]);
caxis([min(z), max(z)]);
% caxis([15, 20]);
zlim([13 23]);

%% Fit Models w/ Subject as Categorical Variable
% Left Hemisphere
load YCS_data2;     load YCS_eloc_MNI;
load YCQ_data_left; load YCQ_eloc_MNI_left;
load YCY_data_left; load YCY_eloc_MNI_left;
load YDD_data;      load YDD_eloc_MNI; 
load YDF_data_left; load YDF_eloc_MNI_left; 
load YDJ_data;      load YDJ_eloc_MNI; 
load YDH_data_left; load YDH_eloc_MNI_left; 

subid1 = 'YCS'; ID_1 = cellstr(repmat(subid1, [size(YCS_data, 1) 1])); 
subid2 = 'YCQ'; ID_2 = cellstr(repmat(subid2, [size(YCQ_data_left, 1) 1])); 
subid3 = 'YCY'; ID_3 = cellstr(repmat(subid3, [size(YCY_data_left, 1) 1])); 
subid4 = 'YDD'; ID_4 = cellstr(repmat(subid4, [size(YDD_data, 1) 1])); 
subid5 = 'YDF'; ID_5 = cellstr(repmat(subid5, [size(YDF_data_left, 1) 1])); 
subid6 = 'YDJ'; ID_6 = cellstr(repmat(subid6, [size(YDJ_data, 1) 1]));
subid7 = 'YDH'; ID_7 = cellstr(repmat(subid7, [size(YDH_data_left, 1) 1]));

% Right Hemisphere
% load YCP_data; load YCP_eloc_MNI; 
% load YCQ_data; load YCQ_eloc_MNI; 
% load YCU_data; load YCU_eloc_MNI; 
% load YCW_data; load YCW_eloc_MNI; 
% load YCY_data; load YCY_eloc_MNI; 
% load YDF_data; load YDF_eloc_MNI; 
% load YDH_data; load YDH_eloc_MNI; 
% 
% subid1 = 'YCP'; ID_1 = cellstr(repmat(subid1, [size(YCP_data, 1) 1])); 
% subid2 = 'YCQ'; ID_2 = cellstr(repmat(subid2, [size(YCQ_data, 1) 1])); 
% subid3 = 'YCU'; ID_3 = cellstr(repmat(subid3, [size(YCU_data, 1) 1])); 
% subid4 = 'YCW'; ID_4 = cellstr(repmat(subid4, [size(YCW_data, 1) 1])); 
% subid5 = 'YCW'; ID_5 = cellstr(repmat(subid5, [size(YCY_data, 1) 1])); 
% subid6 = 'YDF'; ID_6 = cellstr(repmat(subid6, [size(YDF_data, 1) 1])); 
% subid7 = 'YDH'; ID_7 = cellstr(repmat(subid7, [size(YDH_data, 1) 1]));

subIDs = [ID_1; ID_2; ID_3; ID_4; ID_5; ID_6; ID_7];

% Create Table of Variables
var_names = {'Subjects', 'Peak_Frequency', 'Peak_Power', 'x_coord', 'y_coord', 'z_coord'};
T = table(subIDs(~rm_idx), peak_frequencies_rm, peak_power_rm, x_coord(~rm_idx), y_coord(~rm_idx), z_coord(~rm_idx), 'VariableNames', var_names);

frequency_mdl_subj = fitlm(T, 'interactions', 'ResponseVar', 'Peak_Frequency',...
                    'PredictorVars',  {'y_coord', 'z_coord', 'Subjects'},...
                    'CategoricalVar', {'Subjects'}, ... 
                    'RobustOpts', 'on'); 
        
power_mdl_subj = fitlm(T, 'interactions', 'ResponseVar', 'Peak_Power',...
                 'PredictorVars',  {'y_coord', 'z_coord', 'Subjects'},...
                 'CategoricalVar', {'Subjects'}, ... 
                 'RobustOpts', 'on');  

[~, frequency_mdl_pval] = corrcoef(peak_frequencies_rm, frequency_mdl_subj.Fitted, 'rows', 'complete');          
[~, power_mdl_pval]     = corrcoef(peak_power_rm, power_mdl_subj.Fitted, 'rows', 'complete');            

%% JODAL - LOAD Data and Find Freq Bands 
load('INSULA_ALL_LH.mat', 'INSULA_ALL_LH'); 
load('INSULA_ALL_LH.mat', 'INSULA_ALL_LH_coords');
load('INSULA_ALL_RH.mat', 'INSULA_ALL_RH'); 
load('INSULA_ALL_RH.mat', 'INSULA_ALL_RH_coords');

raw_data = INSULA_ALL_LH;
eloc = INSULA_ALL_LH_coords;
x_coord  = eloc(:,1); 
y_coord  = eloc(:,2); 
z_coord  = eloc(:,3); 

% parameters
srate = 2000;
params.srate = srate; % sampling rate
params.wavefreqs = 1:0.1:50; % frequency
params.wavecycles = 6; % wavelet cycles
params.local_winsize_sec = 10; % window size (secs)

% compute psd and fit
[~, pow] = multiphasevec2(params.wavefreqs, mean(raw_data), params.srate, params.wavecycles);
logPow   = log(mean(pow, 2)); 
logFreqs = log(params.wavefreqs);

[b, output] = robustfit(logFreqs, logPow');
[pks, locs] = findpeaks(output.resid);
peak_freqs = params.wavefreqs(locs)'; % select peak frequencies
fitLine = b(1) + b(2) .* logFreqs;

% compute bands (3 Hz on each size)
freq_radius = 3;
lowfreq_lim = peak_freqs - freq_radius;
highfreq_lim = peak_freqs + freq_radius;
bands = [lowfreq_lim, highfreq_lim];
nbands = size(bands,1);

%
figure; 
psd_sem  = std(logPow, 1)/sqrt(size(raw_data, 1));
h = boundedline(double(params.wavefreqs'), double(logPow), double(psd_sem), 'cmap', 'b', 'alpha');
h.LineWidth = 3.0; 
hold on; 
plot(params.wavefreqs, fitLine); 

for band = 2:nbands
    xpoints = [bands(band, 1), bands(band, 1), bands(band, 2), bands(band, 2)];
    ypoints = [0, max(logPow), max(logPow), 0];
    hold on; 
    a = fill(xpoints, ypoints, 'k');
    a.FaceAlpha = 0.05;
end        
set(gcf, 'color', [1 1 1]);
set(gca, 'color', [1 1 1]);
set(gca, 'fontsize', 35);
xlim([1 30]);

%% MULTI OSCILLATORS - PLOT HISTOGRAM
% load('MULTI_OSCILLATORS_ThetaBeta.mat');
figure; 
subplot(2,2,1);
histogram(INSULA_ALL_LH_theta_peakfreqs, 'FaceColor', 'b', 'FaceAlpha', 1.0);
subplot(2,2,2);
histogram(INSULA_ALL_LH_beta_peakfreqs, 'FaceColor', 'b', 'FaceAlpha', 1.0);
subplot(2,2,3);
histogram(INSULA_ALL_RH_theta_peakfreqs, 'FaceColor', 'r', 'FaceAlpha', 1.0);
subplot(2,2,4);
histogram(INSULA_ALL_RH_beta_peakfreqs, 'FaceColor', 'r', 'FaceAlpha', 1.0);

%% MULTI OSCILLATORS - Run MODAL on All Channels
input_data = INSULA_ALL_LH; % specify input data

% parameters
freq_range = [1 50];
params.srate = srate; % sampling rate
params.wavefreqs = freq_range(1):0.1:freq_range(2); % frequency
params.wavecycles = 6; % wavelet cycles
params.local_winsize_sec = 10; % window size (secs)

% operation
nchannels = size(input_data, 1); 
clear multi_oscillators n_oscillators;
tic
for chan = 1:nchannels
    disp(['Computing MODAL for channel ', num2str(chan), ' ...']);
    [~, bands] = MODAL(input_data(chan, :), params); 
    multi_oscillators(chan).bands = bands; % store selected frequency bands
    multi_oscillators(chan).center_freq = mean(bands);
    n_oscillators(chan,:) = size(multi_oscillators(chan).bands, 1); % compute the number of 'oscillators'
end
toc; 

%% MULTI OSCILLATORS - Load and Plot
load('MULTI_OSCILLATORS.mat');
figure; 
subplot(1,2,1);
histogram(n_oscillators_LH, 'FaceColor', 'r', 'FaceAlpha', 1.0);
subplot(1,2,2);
histogram(n_oscillators_RH, 'FaceColor', 'b', 'FaceAlpha', 1.0);

%% %% m.MODAL (revised) - Find Peaks Across Specified Bands
clear peak_frequencies peak_power all_power;
% load('YDD_data'); % load pre-processed data
raw_data = INSULA_ALL_RH;

% Compute MODAL first
% [frequency_sliding, bands, bandpow, bandphases, pow] = MODAL(mean(INSULA_ALL_RH), params); 

srate = 2000; 
% freq_range = bands(1,:); 
freq_range = [1 50]; 
nchannels = size(raw_data, 1); 
params.srate = srate; % sampling rate
params.wavefreqs = freq_range(1):0.1:freq_range(2); % frequency
params.wavecycles = 6; % wavelet cycles
params.local_winsize_sec = 10; % window size (secs)
color_set = distinguishable_colors(nchannels,'k'); 

% compute PSD
figure; 
for chan = 1:nchannels
    disp(['Computing PSD from ', num2str(freq_range(1)), '-', num2str(freq_range(2)), ' Hz for channel ', num2str(chan), ' ...']);
    [phase, pow] = multiphasevec2(params.wavefreqs, raw_data(chan, :), params.srate, params.wavecycles);
    logPow   = log(mean(pow,2)); % vector of power values from wavelet averaged across time for each frequency 
    logFreqs = log(params.wavefreqs);
    [b, output] = robustfit(logFreqs, logPow');
    fit_line = b(1) + b(2) .* logFreqs; % 1/f
    [pks, locs] = findpeaks(output.resid);
    peak_freqs = params.wavefreqs(locs)'; % select peak frequencies
    [power, max_power_idx] = max(pks); % find index of maximum power
    peak_frequency         = peak_freqs(max_power_idx); % peak_frequency = frequency with greatest power    
    power_peak             = logPow(max_power_idx); % power level of peak frequency    
    
    if  ~isempty(peak_frequency)
        peak_frequencies(:, chan) = peak_frequency; % stores only peak frequencies
        peak_power(:, chan) = power_peak;
    else 
        peak_frequencies(:, chan) = NaN; % use zero if there was no peak
        peak_power(:, chan) = NaN;
    end
    all_power(chan, :) = logPow; 
    plot(params.wavefreqs, logPow, 'Color', color_set(chan,:), 'LineWidth', 5.0); 
%     hold on; 
%     plot(params.wavefreqs, fitLine); 
end

% peak_freqs             = freqs(locs)'; % select peak frequencies
% [power, max_power_idx] = max(pks); % find index of maximum power
% peak_frequency         = peak_freqs(max_power_idx); % peak_frequency = frequency with greatest power
% peak_frequencies(:,p)  = peak_frequency; % stores only peak frequencies

%%