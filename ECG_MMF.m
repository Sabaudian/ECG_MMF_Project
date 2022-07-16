%% Load Dataset and apply butterworth filter on it
clc; % Clear the command window.
close all; % Close all figures (except those of imtool.)
clear; % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 15;
markerSize = 8;
% load all the data 
G = 200; % Gain
Fs = 360; % [Hz]
L = 3600; % lenght of ECG signals
T = linspace(0,L/Fs,L); % time axis
F = linspace(-Fs/2, Fs/2, L); % Frequency axis
files = dir(fullfile("dataset/","*.mat")); % all dataset files
numData = numel(files); % number of data
ECGs = zeros(numData,L); % prealloc
% load and store data
for i = 1:numData
    load(fullfile("dataset/",files(i).name)); % load all data
    ECGs(i,:) = val/G;
end
% Plot the signal/s you want
figure(1); plot(T, ECGs(1,:)); grid on;
title("ECG Signal","FontSize",fontSize); 
xlabel("Time (sec)", "FontSize", fontSize); 
ylabel("voltage [mV]", "FontSize", fontSize);
% Define a Butterworth Filter
[b,a] = butter(3,[1 30]/(Fs/2),"bandpass"); 
FLT_ECGs = zeros(numData,L); % prealloc 
% filter all signals
for i = 1:numData
    FLT_ECGs(i,:) = filtfilt(b,a,ECGs(i,:)); 
end
% Plot the Filterd signals you want
figure(2); plot(T, FLT_ECGs(1,:)); grid on;
title("Clean ECG signal","FontSize",fontSize); 
xlabel("Time (sec)", "FontSize", fontSize); 
ylabel("voltage [mV]", "FontSize", fontSize);
clear a b files val;

%% Add random generated Noise to all the signals 
close all;
% Prealloc 
NS_ECGs = zeros(numData, L); 
SNR = zeros(1,numData); 
% Add noise to all signals
for i = 1:numData
    % generate random snr for all signals
    SNR(i) = 1 +  (10 - 1).*rand(1);
    % additive white gaussian noise
    NS_ECGs(i,:) = awgn(FLT_ECGs(i,:),SNR(i),'measured');
end
% Plot the Noisy Signal of the signals
figure(1); plot(T, NS_ECGs(1,:), "b-"); grid on;
title("Noisy ECG Signal", "FontSize", fontSize);  
xlabel("Time (sec)", "FontSize", fontSize); 
ylabel("Voltage (Hz)", "FontSize", fontSize);
clear SNR;

%% Add random generated Baseline drift to all signals
close all;
% prealloc
DFT_ECGs = zeros(numData,L);
% Call funtion to generate drift
drift = GenDrift(numData,L);
% Apply drift to all signals
for i = 1:numData
    DFT_ECGs(i,:) = NS_ECGs(i,:) + drift(i,:);
end
% Plot the Filterd signal/s you want
figure(1); plot(T, DFT_ECGs(1,:), "b-"); grid on;
title("Drifted ECG signal", "FontSize", fontSize); 
xlabel("Time (sec)", "FontSize", fontSize);  
ylabel("Voltage (Hz)", "FontSize", fontSize);
clear drift;

%% Baseline correction
close all;
% Defining the two structuring element
Bo = ones(1,0.2*Fs+1); % Bo = strel('line',0.2*Fs,0);
Bc = ones(1,round(1.5*0.2*Fs+1)); % Bc = strel('line',1.5*0.2*Fs,0);
% Prealloc
peaksSuppression = zeros(numData,L);
pitsSuppression = zeros(numData,L);
detectedDrift = zeros(numData,L);
Correction = zeros(numData,L);
finalBaseline = zeros(numData,L);
% Opening and Closing application to all signals
for i=1:numData
    % Opening: erosion B dilatation B
    peaksSuppression(i,:) = opening(DFT_ECGs(i,:), Bo); % peaksSuppression(i,:) = imopen(DFT_ECGs(i,:), Bo);
    % closing: dilatation B erosion B
    pitsSuppression(i,:) = closing(DFT_ECGs(i,:), Bc); % pitsSuppression(i,:) = imclose(DFT_ECGs(i,:), Bo);
end
% Plot the representation of Op. and Cls. operation 
figure(1); subplot(2,1,1); hold on;
plot(T, peaksSuppression(1,:),"g-","LineWidth",3);
plot(T, DFT_ECGs(1,:), "b-","LineWidth",0.5);
title("Opening operation", "FontSize", fontSize); 
xlabel("Time (sec)", "FontSize", fontSize);  
ylabel("Voltage (Hz)", "FontSize", fontSize);
legend("Opening","Signal");
grid on; hold off;
subplot(2,1,2); hold on;
plot(T, pitsSuppression(1,:),"g-","LineWidth",3);
plot(T, DFT_ECGs(1,:), "b-","LineWidth",0.5);
title("Closing operation", "FontSize", fontSize); 
xlabel("Time (sec)", "FontSize", fontSize);  
ylabel("Voltage (Hz)", "FontSize", fontSize);
legend("Closing","Signal");
grid on; hold off;
% Detection and Correction of the Wandering Baseline 
% Apply Op. and Cls. to detect the drift
for i=1:numData
    peaksSuppression(i,:) = opening(DFT_ECGs(i,:), Bo); % peaksSuppression(i,:) = imopen(DFT_ECGs(i,:), Bo);
    pitsSuppression(i,:) = closing(peaksSuppression(i,:), Bc); % pitsSuppression(i,:) = imclose(peaksSuppression(i,:), Bc);
    % Detected drift of all sinals
    detectedDrift(i,:) = pitsSuppression(i,:);
end
% Correction subtracting the drift from signals
for i=1:numData
    % Signal With Baseline drift corrected
    Correction(i,:) = DFT_ECGs(i,:) - detectedDrift(i,:);
    % Finale Baseline result
    finalBaseline(i,:) = closing(opening(Correction(i,:),Bo),Bc); % finalBaseline(i,:) = imclose(imopen(Correction(i,:), Bo),Bc);
end
% Plot Corrected Signal and Baseline
figure(2); subplot(2,1,1); hold on;
plot(T, detectedDrift(1,:),"g-","LineWidth",3);
plot(T, DFT_ECGs(1,:), "b-","LineWidth",0.5);
title("Detected Baseline Drift", "FontSize", fontSize); 
xlabel("Time (sec)", "FontSize", fontSize);  
ylabel("Voltage (Hz)", "FontSize", fontSize);
legend("Baseline","Signal");
grid on; hold off;
subplot(2,1,2); hold on;
plot(T, finalBaseline(1,:),"g-","LineWidth",3);
plot(T, Correction(1,:), "b-","LineWidth",0.5);
title("Drift Correction", "FontSize", fontSize); 
xlabel("Time (sec)", "FontSize", fontSize);  
ylabel("Voltage (Hz)", "FontSize", fontSize);
legend("Baseline","Signal");
grid on; hold off;

%% Noise Suppression with MMF Algorithm
close all;
% Prealloc
dilatation_1 = zeros(numData,L);
erosion_1 = zeros(numData,L);
dilatation_2 = zeros(numData,L);
erosion_2 = zeros(numData,L);
MMF_denoise = zeros(numData,L);
% Defining the two structuring element
B1 = [0 1 5 1 0]; B2 = [1 1 1 1 1];
% Apply the Algorithm for all signals
% 1/2*((Fbc dilatation B1 erosion B2) + (Fbc erosion B1 dilatation B2))
for i = 1:numData
    dilatation_1(i,:) = dilatation(Correction(i,:), B1); % dilatation
    erosion_1(i,:) = erosion(dilatation_1(i,:), B2); % erosion
    erosion_2(i,:) = erosion(Correction(i,:), B1); % erosion
    dilatation_2(i,:) = dilatation(erosion_2(i,:), B2); % dilatation
    MMF_denoise(i,:) = (erosion_1(i,:) + dilatation_2(i,:))/2; % average
end
% Plot denoised signal with MMF algorithm
figure(1); plot(T, MMF_denoise(1,:), "b-"); grid on;
title("MMF Denoised Ecg signal", "FontSize", fontSize);  
xlabel("Time (sec)", "FontSize", fontSize); 
ylabel("Voltage (Hz)", "FontSize", fontSize);
clear dilatation_1 dilatation_2 erosion_1 erosion_2;

%% Noise Suppressione with MF Algorithm
% Prealloc
MF_op_cl = zeros(numData,L);
MF_cl_op = zeros(numData,L);
MF_denoise = zeros(numData,L);
% Defining the structuring element
B = [0 1 5 1 0];
% Apply the Algorithm for all signals
% 1/2*((Fbc opening B closing B) + (Fbc closing B opening B))
for i = 1:numData
    % op. and cl.
    MF_op_cl(i,:) = closing(opening(Correction(i,:),B),B); 
    % cl. and op.
    MF_cl_op(i,:) = opening(closing(Correction(i,:),B),B); 
    % average
    MF_denoise(i,:) = (MF_op_cl(i,:) + MF_cl_op(i,:))/2; 
end
% Plot denoised signal with Chu's MF algorithm
figure(2); plot(T, MF_denoise(1,:), "b-"); grid on;
title("CHU's MF Denoised Ecg signal", "FontSize", fontSize);  
xlabel("Time (sec)", "FontSize", fontSize); 
ylabel("Voltage (Hz)", "FontSize", fontSize);
clear MF_op_cl MF_cl_op;

%% Evaluation of the algorithms
close all;
% Prealloc
MMF_NSR = zeros(1,numData);
MF_NSR = zeros(1,numData);
MMF_SDR = zeros(1,numData);
MF_SDR = zeros(1,numData);
% Compute NSR and SDR of the two method
for i = 1:numData
    % NSR (Noise Suppression Ratio) bigger the better
    MMF_NSR(i) = sum(abs(fft(MMF_denoise(i,:))))/sum(abs(fft(FLT_ECGs(i,:)))); % MMF's NSR of all signals
    MF_NSR(i) = sum(abs(fft(MF_denoise(i,:))))/sum(abs(fft(FLT_ECGs(i,:)))); % MF's NSR of all signals
    % SDR (Signal Distortion Ratio) smaller the better
    MMF_SDR(i) = sum(abs(fft(FLT_ECGs(i,:) - MMF_denoise(i,:))))/sum(abs(fft(MMF_denoise(i,:)))); % MMF's SDR of all signal
    MF_SDR(i) = sum(abs(fft(FLT_ECGs(i,:) - MF_denoise(i,:))))/sum(abs(fft(MF_denoise(i,:)))); % MF's SDR of all signal
end
% Graphical Representation of the NSR
figure(1); hold on; grid on;
plot(MF_NSR,'-^r','LineWidth',1,'MarkerSize',markerSize);
plot(MMF_NSR,'-ob','LineWidth',1,'MarkerSize',markerSize);
hold off; xticks(0:50);
title("Comparison of NSRs", "FontSize", fontSize);
xlabel("Signals", "FontSize", fontSize); 
ylabel("NSR", "FontSize", fontSize);
legend("Chu's MF Algorithm", "MMF Algorithm", "FontSize", fontSize);

% Graphical Representation of the SDR
figure(2); hold on; grid on;
plot(MF_SDR,'-^r','LineWidth',1,'MarkerSize',markerSize);
plot(MMF_SDR,'-ob','LineWidth',1,'MarkerSize',markerSize);
hold off; xticks(0:50);
title("comparison of SDRs", "FontSize", fontSize);
xlabel("Signals", "FontSize", fontSize); 
ylabel("SDR", "FontSize", fontSize);
legend("Chu's MF Algorithm","MMF Algorithm", "FontSize", fontSize);

%% Evaluation changing the dimension of the structuring element
close all;
% Define the dim. of strel
N = 100;
% Prealloc
S1 = cell(1,N); % strel one
S2 = cell(1,N); % strel two
dataset = cell(1,N); % dataset as cell array
snr = zeros(1,N);
% Generate Structuring Elements
for i = 1:N
    S1(i) = {GenStrel(i)}; % strel 1
    S2(i) = {ones(1,length(S1{i}))}; % strel 2
    dataset(i) = {FLT_ECGs};
    snr(i) = 1 + (30 - 1).*rand(1); % define random snr
end
% Generate increasing noise
snr = sort(snr,"descend");
noisy_dataset = cell(1,N);
% Add noise to dataset
for i = 1:N
    noisy_dataset(i) = {awgn(dataset{i},snr(i),'measured')};
end

% Apply MMF algorithm with various strel length
curr_mmf = zeros(numData,L);
all_mmf = cell(1,N);
dl_er = zeros(numData,L);
er_dl = zeros(numData,L);
for j=1:N
    s1 = S1{j};
    s2 = S2{j};
    for i=1:numData  % Compute the Algorithm
        dl_er(i,:) = erosion(dilatation(noisy_dataset{j}(i,:), s1), s2);
        er_dl(i,:) = dilatation(erosion(noisy_dataset{j}(i,:), s1), s2);
        curr_mmf(i,:) = (dl_er(i,:) + er_dl(i,:))/2;
        all_mmf(j) = {curr_mmf}; 
    end
end
clear curr_mmf dl_er er_dl; % clear no more necessary variable

% Prealloc
curr_nsr = zeros(1,numData);
all_mmf_nsr = cell(1,N);
curr_sdr = zeros(1,numData);
all_mmf_sdr = cell(1,N);
% Compute NSRs and SDRs
for j=1:N
    current_dataset = all_mmf{j}; % select current dataset
    for i=1:numData
        % MMF NSR:
        curr_nsr(i) = sum(abs(fft(current_dataset(i,:))))/sum(abs(fft(FLT_ECGs(i,:)))); % current dataset NSR values
        all_mmf_nsr(j) = {curr_nsr}; % all values of NSRs per dataset defined as strel grows
        % MMF SDR:
        curr_sdr(i) = sum(abs(fft(FLT_ECGs(i,:)-current_dataset(i,:))))/sum(abs(fft(current_dataset(i,:))));
        all_mmf_sdr(j) = {curr_sdr}; % nsr di tutti e 48 i messaggi, uno per ogni strel
    end
end
clear current_dataset curr_nsr curr_sdr; % clear no more necessary variable

% Prealloc
op_cl = zeros(numData,L);
cl_op = zeros(numData,L);
curr_mf = zeros(numData,L);
all_mf = cell(1,N);
% Apply Chu's MF algorithm with various strel length
for j=1:N
    s1 = S1{j};
    % Compute the MF Algorithm
    for i=1:numData  
        op_cl(i,:) = closing(opening(noisy_dataset{j}(i,:),s1),s1);
        cl_op(i,:) = opening(closing(noisy_dataset{j}(i,:),s1),s1);
        curr_mf(i,:) = (op_cl(i,:) + cl_op(i,:))/2;
        all_mf(j) = {curr_mf}; 
    end
end
clear curr_mf op_cl cl_op all_mmf s1 s2; % clear no more necessary variable

% Prealloc
curr_nsr = zeros(1,numData);
all_mf_nsr = cell(1,N);
curr_sdr = zeros(1,numData);
all_mf_sdr = cell(1,N);
% Compute NSRs and SDRs
for j=1:N
    current_dataset = all_mf{j}; % select current dataset
    for i=1:numData
        % MMF NSR:
        curr_nsr(i) = sum(abs(fft(current_dataset(i,:))))/sum(abs(fft(FLT_ECGs(i,:)))); % current dataset NSR values
        all_mf_nsr(j) = {curr_nsr}; % all values of NSRs per dataset defined as strel grows
        % MMF SDR:
        curr_sdr(i) = sum(abs(fft(FLT_ECGs(i,:)-current_dataset(i,:))))/sum(abs(fft(current_dataset(i,:))));
        all_mf_sdr(j) = {curr_sdr}; % nsr di tutti e 48 i messaggi, uno per ogni strel
    end
end
clear current_dataset curr_nsr curr_sdr all_mf; % clear no more necessary variable

% Defining mean of all mmf and mf result and confront it 
mean_mmf_nsr = cellfun(@(x) mean(x, "all"), all_mmf_nsr);
mean_mmf_sdr = cellfun(@(x) mean(x, "all"), all_mmf_sdr);
mean_mf_nsr = cellfun(@(x) mean(x, "all"), all_mf_nsr);
mean_mf_sdr = cellfun(@(x) mean(x, "all"), all_mf_sdr);
% Evaluation
% Graphical Representation of the NSR
figure(1); hold on; grid on;
plot(mean_mf_nsr,'-^r','LineWidth',1);
plot(mean_mmf_nsr,'-ob','LineWidth',1);
hold off;
get(gca,'XTickLabel'); set(gca,'XTickLabel',(1:2:2*N-1));
xticks(1:2*N-1);
title("NSRs as Strel grows", "FontSize", fontSize);
xlabel("strel growing", "FontSize", fontSize); 
ylabel("NSR", "FontSize", fontSize);
legend("Chu's MF Algorithm", "MMF Algorithm");
% Graphical Representation of the SDR
figure(2); hold on; grid on;
plot(mean_mf_sdr,'-^r','LineWidth',1);
plot(mean_mmf_sdr,'-ob','LineWidth',1);
get(gca,'XTickLabel'); set(gca,'XTickLabel',(1:2:2*N-1));
hold off; xticks(1:2*N-1);
title("SDRs as Strel grows", "FontSize", fontSize);
xlabel("strel growing", "FontSize", fontSize); 
ylabel("SDR", "FontSize", fontSize);
legend("Chu's MF Algorithm","MMF Algorithm");
