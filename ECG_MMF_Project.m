%% First Part: 
% -------------------------------------------------------------- %
% ECG signal simulated by ECGwaveGen function from PhysioNet.org %
% -------------------------------------------------------------- %
bpm = 72; 
sample_frquency = 200; % Hz
time_interval = 10; % seconds
amplitude = 1000; % mV

sim_ecg = ECGwaveGen(bpm,time_interval, sample_frquency, amplitude);

l = length(sim_ecg);
i = (0 : l - 1) / sample_frquency;
plot(i, sim_ecg);
title('Simulated ECG signal');
xlabel ('time [sec]');
ylabel ('voltage [mV]');
grid on

%% Add Noise to the Simulated ECG signal
noisy_sig = awgn(sim_ecg, 20, 'measured');
plot(i, noisy_sig, 'b');
grid on
xlabel('time [sec]');
ylabel('voltage [mV]');
title('Simulated ECG signal with added noise');

%% Add Baseline drift to the Simulated ECG signal
x = linspace(0, 2*pi, l);
A = 0.8;
N = 60;
cos_drift = A * cos(x ./ N);
minDriftOffset = 0; 
maxDriftOffset = 120;
add_drift = cos_drift + linspace(minDriftOffset, maxDriftOffset, length(noisy_sig));
subplot(2,1,1)
plot(i, noisy_sig);
xlabel('time [sec]');
ylabel('voltage [mV]');
title("Noisy ECG Signal");
grid on
grid minor
subplot(2,1,2)
dirty_sig = noisy_sig + add_drift; % Corrupted ECG Signal %
plot(i, dirty_sig)
xlabel('time [sec]');
ylabel('voltage [mV]');
title("Corrupted Simulated ECG Signal (Noise + Baseline Drift)");
grid on
grid minor

%% Baseline Correction on Simulated ECG
% Bo = ones(1, 0.2*fs); % = 0.2*Fs (fs = 360 Hz) %
% Bc = ones(1, 1.5*0.2*fs); % = 1.5*o %

bo = strel('line', 0.2*sample_frquency, 0);
bc = strel('line', 1.5*0.2*sample_frquency, 0);

op = imopen(dirty_sig, bo); % suppress peaks %
figure(1)
subplot(2,1,1)
hold on
plot(i, op, 'g', 'LineWidth', 3)
plot(i, dirty_sig, 'b', 'LineWidth', 0.5)
xlabel('time [sec]');
ylabel('voltatge [mV]');
title('Opening operation');
legend('Opening','ECG Signal');
grid on
hold off

cls = imclose(dirty_sig, bc); % suppress pits %
%figure(2)
subplot(2,1,2)
hold on
plot(i, cls, 'g', 'LineWidth', 3);
plot(i, dirty_sig, 'b', 'LineWidth', 0.5)
xlabel('time [sec]');
ylabel('voltage [mV]');
title('Closing operation');
legend('Closing','ECG Signal');
grid on
hold off

detect_drift = imclose(op, bc); % Baseline drift detected %
figure(2)
subplot(2,1,1)
hold on
plot(i, detect_drift, 'g', 'LineWidth', 3)
plot(i, dirty_sig, 'b', 'LineWidth', 0.5)
xlabel('time [sec]');
ylabel('voltage [mV]');
title('Baseline Drift detected in the ECG');
legend('Baseline drift', 'ECG signal');
grid on
hold off

drift_correction = dirty_sig - detect_drift; % drift deleted %
subplot(2,1,2)
hold on
final_op = imopen(drift_correction, bo);
final_baseline = imclose(final_op, bc);
plot(i, final_baseline, 'g', 'LineWidth', 3)
plot(i, drift_correction, 'b', 'LineWidth', 0.5)
grid on
title('Baseline Corrected ECG Signal');
xlabel('time [sec]');
ylabel('voltage [mV]');
title("ECG signal after Baseline drift correction");
legend('Corrected Baseline', 'ECG signal');
hold off

%% MMF Algorithm: Noise Suppression 
% f = 1/2 * (f_bc • B_pair + f_bc o B_pair)
% f_bc: signal after baseline correction
% B_pair = {B1, B2}
b1 = [0 1 5 1 0]; % triangular shape
b2 = [1 1 1 1 1]; % linear shape [1 1 1 1 1]

% first addend: f_bc dilatation B1 erosion B2 %
v_1 = imdilate(drift_correction, b1); % "expansion" %
first_value = imerode(v_1, b2); % "shrinking" %

% Second addend: f_bc erosion B1 dilataion B2
v_2 = imerode(drift_correction, b1); % "shrinking" %
second_value = imdilate(v_2, b2); % "expansion" %

denoise_signal = (first_value + second_value) / 2; 

subplot(2,1,1)
plot(i, denoise_signal, 'b');
title("ECG signal conditioning with MMF algorithm");
xlabel("time [sec]");
ylabel("voltage [mV]");
grid on
grid minor

%% CHU's MF Algorithm: Noise Suppression on Simulated ECG
% Impulsive noise suppression is performed
% by processing the data through a sequence
% of opening and closing operations. 
% The result from this step is the average 
% of the two estimates. 
% Only one structuring element (B) is used.

% f = 1/2 * ((f_in o B • B) + (f_in • B o B))

b = [0 1 5 1 0]; % triangular shape

% first part: opening and closing 
mf_op_1 = imopen(drift_correction, b);
mf_cl_1 = imclose(mf_op_1, b);

% second part: closing and opening
mf_cl_2 = imclose(drift_correction, b);
mf_op_2 = imopen(mf_cl_2, b);

% average
mf_noise_suppression = (mf_op_1 + mf_cl_2) / 2;

subplot(2,1,2)
plot(i, mf_noise_suppression);
title("ECG signal conditioning with CHU's MF Algorithm");
xlabel("time [sec]");
ylabel("voltage [mV]");
grid on
grid minor

%% Evalutation on Simulated ECG signal
% BCR (Baseline Correction Ratio) bigger the better
bcr_mmf = sum(abs(drift_correction)) / sum(abs(detect_drift));
fprintf('> BCR = %f\n\n', bcr_mmf);

% NSR (Noise Suppression Ratio) bigger the better
nsr_mmf = sum(abs(denoise_signal)) / sum(abs(sim_ecg));
fprintf('> NSR_MMF = %f | ', nsr_mmf);

nsr_mf = sum(abs(mf_noise_suppression)) / sum(abs(sim_ecg));
fprintf('NSR_MF = %f\n\n', nsr_mf);

% SDR (Signal Distortion Ratio) smaller the better
sdr_mmf = sum(abs(sim_ecg - denoise_signal)) / sum(abs(denoise_signal));
fprintf('> SDR_MMF = %f | ', sdr_mmf);

sdr_mf = sum(abs(sim_ecg - mf_noise_suppression)) / sum(abs(mf_noise_suppression));
fprintf('SDR_MF = %f\n\n', sdr_mf);

%% Second Part:
% -------------------------------------------------------------- %
% Real ECG signal recovered from PhysioNet.org database
% Load and plot the original ECG signal
% -------------------------------------------------------------- %
load('100m.mat');
original_signal = val(1,:);
fs = 360; % Hz
% T = 1/fs; 
L = length(original_signal);
t = (0 : L - 1) / fs;

plot(t, original_signal);
title('plot of the original ECG signal');
xlabel ('time [sec]');
ylabel ('voltage [mV]');
grid on

%% Baseline Correction 
Bo = strel('line', 0.2*fs, 0);
Bc = strel('line', 1.5*0.2*fs, 0);

opening = imopen(original_signal, Bo); % suppress peaks %
figure(1)
subplot(2,1,1)
hold on
plot(t, opening, 'g', 'LineWidth', 3)
plot(t, original_signal, 'b', 'LineWidth', 0.5)
xlabel('time [sec]');
ylabel('voltatge [mV]');
title('Opening operation');
legend('Opening','ECG Signal');
grid on
hold off

closing = imclose(original_signal, Bc); % suppress pits %
%figure(2)
subplot(2,1,2)
hold on
plot(t, closing, 'g', 'LineWidth', 3);
plot(t, original_signal, 'b', 'LineWidth', 0.5)
xlabel('time [sec]');
ylabel('voltage [mV]');
title('Closing operation');
legend('Closing','ECG Signal');
grid on
hold off

correction = imclose(opening, Bc); % Baseline drift detected %
figure(2)
subplot(2,1,1)
hold on
plot(t, correction, 'g', 'LineWidth', 3)
plot(t, original_signal, 'b', 'LineWidth', 0.5)
xlabel('time [sec]');
ylabel('voltage [mV]');
title('Baseline Drift detected in the ECG');
legend('Baseline drift', 'ECG signal');
grid on
hold off

BaselineDrift_Correction = original_signal - correction; % drift deleted %
subplot(2,1,2)
hold on
check_op = imopen(BaselineDrift_Correction, Bo);
check_baseline = imclose(check_op, Bc);
plot(t, check_baseline, 'g', 'LineWidth', 3)
plot(t, BaselineDrift_Correction, 'b', 'LineWidth', 0.5)
grid on
title('Baseline Corrected ECG Signal');
xlabel('time [sec]');
ylabel('voltage [mV]');
title("ECG signal after Baseline drift correction");
legend('Corrected Baseline', 'ECG signal');
hold off

%% MMF Algorithm: Noise Suppression 
B1 = [0 1 5 1 0]; % triangular shape
B2 = [1 1 1 1 1]; % linear shape [1 1 1 1 1]

% first addend: f_bc dilatation B1 erosion B2 %
value_1 = imdilate(BaselineDrift_Correction, B1); % "expansion" %
first_addend = imerode(value_1, B2); % "shrinking" %

% Second addend: f_bc erosion B1 dilataion B2
value_2 = imerode(BaselineDrift_Correction, B1); % "shrinking" %
second_addend = imdilate(value_2, B2); % "expansion" %

Denoise_signal = (first_addend + second_addend) / 2; 

subplot(2,1,1)
plot(t, Denoise_signal, 'b');
title("ECG signal conditioning with MMF algorithm");
xlabel("time [sec]");
ylabel("voltage [mV]");
grid on
grid minor

%% CHU's MF Algorithm: Noise Suppression 
B = [0 1 5 1 0]; % triangular shape

% first part: opening and closing 
mf_opening_1 = imopen(BaselineDrift_Correction, B);
mf_closing_1 = imclose(mf_opening_1, B);

% second part: closing and opening
mf_closing_2 = imclose(BaselineDrift_Correction, B);
mf_opening_2 = imopen(mf_closing_2, B);

% average
mf_denoise = (mf_opening_1 + mf_closing_2) / 2;

subplot(2,1,2)
plot(t, mf_denoise);
title("ECG signal conditioning with CHU's MF Algorithm");
xlabel("time [sec]");
ylabel("voltage [mV]");
grid on
grid minor

%% Evalutation
% BCR (Baseline Correction Ratio) bigger the better
BCR_MMF = sum(abs(correction)) / sum(abs(BaselineDrift_Correction));
fprintf('> BCR = %f\n\n', BCR_MMF);

% NSR (Noise Suppression Ratio) bigger the better
NSR_MMF = sum(abs(Denoise_signal)) / sum(abs(original_signal));
fprintf('> NSR_MMF = %f | ', NSR_MMF);

NSR_MF = sum(abs(mf_denoise)) / sum(abs(original_signal));
fprintf('NSR_MF = %f\n\n', NSR_MF);

% SDR (Signal Distortion Ratio) smaller the better
SDR_MMF = sum(abs(original_signal - Denoise_signal)) / sum(abs(Denoise_signal));
fprintf('> SDR_MMF = %f | ', SDR_MMF);

SDR_MF = sum(abs(original_signal - mf_denoise)) / sum(abs(mf_denoise));
fprintf('SDR_MF = %f\n\n', SDR_MF);
