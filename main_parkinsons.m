%parse data
%use init_processing along with parseData function

%% main_parkinsons
%Nx7 matrices, saved by ID
all_subjects = ["001A", "002A","004A", "010A", "115A", "118A", "120A", "215A", "218A", "220A", "031B", "079B", "111B", "211B", "121B", "221B"]';
totalsubjects = length(all_subjects);
peaknumbers = [];

%%
for subject = 1:totalsubjects
    load(strcat('kav',extractBefore(all_subjects(subject),4),'_acc.mat'));
    load(strcat('kav',extractBefore(all_subjects(subject),4),'_gyro.mat'));
    load(strcat('kav',extractBefore(all_subjects(subject),4),'_orien.mat'));
    Tentries_acc = table2array(Tentries_acc); Tentries_gyro = table2array(Tentries_gyro); Tentries_orien = table2array(Tentries_orien);
    matrix = interpolate_data(Tentries_acc(1:end-1, :), Tentries_gyro(1:end-1, :), Tentries_orien(1:end-1, :));
    save(strcat('kav', all_subjects(subject), '_main.mat'), 'matrix');
end
%figure;

%%
clf;
for subject = 1:totalsubjects
    id = char(all_subjects(subject));
    load(strcat('kav',id,'_main.mat'));
    
%% apply low pass filter to smooth data
    sfq = 100; %sampling frequency in Hz
    cfq =10; %cutoff frequency in Hz
    low_cutoff = cfq/(sfq/2); %high cutoff -> more inclusive 
    [b,a] = butter(1,low_cutoff, 'low'); %filter params for cutoff of 0.2
    data_acc_sm = zeros(size(matrix));
    data_acc_sm(:,2:end) = filter(b,a,matrix(: ,2:end)); %applying filter to all accel data
    %plot parkinsons vs non parkinsons 
    if id(4) == 'A'
        figure(1); set(gcf, 'name', 'PD Raw and Filtered Signals');
        subplot(2, 5, subject);
        plot(matrix(:, 1),matrix(:, 2),'b',matrix(:,1),data_acc_sm(:, 2),'r');        
        title(strcat('kav',id));
        xlabel('time (ms)'); ylabel('acceleration(m/s^2)');
    else
        figure(2); set(gcf, 'name', 'No-PD Raw and Filtered Signals');
        subplot(2, 3, subject-10);
        plot(matrix(:, 1),matrix(:, 2),'b',matrix(:,1),data_acc_sm(:, 2),'r');        
        title(strcat('kav',id));
        xlabel('time (ms)'); ylabel('acceleration(m/s^2)');
    end
    matrix(:, 2:end) = data_acc_sm(:, 2:end);
%% spectrograms
    %plot spectrograms
    frequencyLimits = [0 100]/pi; %Normalized frequency (*pi rad/sample)
    leakage = 0.2;
    overlapPercent = 50;
    
    if id(4) == 'A'
        figure(3); set(gcf, 'name', 'PD Spectrogram'); 
        subplot(2, 5, subject);
        pspectrum(data_acc_sm(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
        'Leakage',leakage, 'OverlapPercent',overlapPercent);
        title(strcat('kav',id));
        xlabel('time (ms)'); ylabel('acceleration(m/s^2)');
    else
        figure(4); set(gcf, 'name', 'No-PD Spectrogram'); 
        subplot(2, 3, subject-10);
        pspectrum(data_acc_sm(:, 2), 'spectrogram','FrequencyLimits',frequencyLimits, ...
        'Leakage',leakage, 'OverlapPercent',overlapPercent);
        title(strcat('kav',id));
        xlabel('time (ms)'); ylabel('acceleration(m/s^2)');
    end

%% fourier transform
 % fourier transform 
    Y = fft(matrix(:, 2));
    L = size(matrix(:, 2), 1);
    %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P2 = abs(Y); P1 = P2(1:L/2+1);
    %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
    f = (1000/15)*[0:(L/2)]/L;
    
    if id(4) == 'A'
        figure(5); set(gcf, 'name', 'PD Single-Sided Amplitude Spectrum of X(t)'); 
        subplot(2, 5, subject);
        plot(f,P1);
        title(strcat('kav',id));
        xlabel('f (Hz)'); ylabel('|P1(f)|');
    else
        figure(6); set(gcf, 'name', 'No-PD Single-Sided Amplitude Spectrum of X(t)'); 
        subplot(2, 3, subject-10);
        plot(f,P1);
        title(strcat('kav',id));
        xlabel('f (Hz)'); ylabel('|P1(f)|');
    end
    
    
%% Energy calculation (e.g. for peak detection)
% energy = data_acc_sm(:,2).^2 + data_acc_sm(:,3).^2 + data_acc_sm(:,4).^2;
%     if id(4) == 'A'
%         figure(7); set(gcf, 'name', 'PD Energy');
%         subplot(2, 5, subject);
%         plot(matrix(:, 1),energy);        
%         title(strcat('kav',id));
%         xlabel('time (ms)'); ylabel('squared acceleration (m/s^2)');
%     else
%         figure(8); set(gcf, 'name', 'No-PD Energy');
%         subplot(2, 3, subject-10);
%         plot(matrix(:, 1),energy);        
%         title(strcat('kav',id));
%         xlabel('time (ms)'); ylabel('squared acceleration (m/s^2)');
%     end
    
%% peak detection
    energy = data_acc_sm(:,2).^2 + data_acc_sm(:,3).^2 + data_acc_sm(:,4).^2;
    [peak, peakLocInds] = findpeaks(energy, 'minPeakHeight', 2.5, 'minPeakDistance', 20);
    peaknumbers = [peaknumbers, length(peak)];
    segmentA = 10;
    segmentB = 0;
    segmentStartIdxs = peakLocInds - segmentA;
    segmentEndIdxs = peakLocInds + segmentB;
    
%     plot(peakLocs,peaks,'r.');
    
    


     [xpeak, xpeakLocInds] = findpeaks(matrix(:,2), 'minPeakHeight', 2.5, 'minPeakDistance', 20);
     segmentA = 10;
     segmentB = 0;
    segmentStartIdxs = xpeakLocInds - segmentA;
    segmentEndIdxs = xpeakLocInds + segmentB;

    if id(4) == 'A'
        figure(9); set(gcf, 'name', 'PD Peak-Detection')
        subplot(2,5,subject);
        plot(matrix(:,1),energy);
        hold on;
        %plot(matrix(:,1),matrix(:,3), matrix(:,1),matrix(:,4));
        timestamps = matrix(:,1);
        peak_timestamps = timestamps(xpeakLocInds);
        plot(peak_timestamps, peak, 'r.');
        plot(xpeakLocs,peaks,'r.');
        title(strcat('kav',id));
        xlabel('time (ms)'); ylabel('acceleration (m/s^2)');
        
        start_seg = [timestamps(segmentStartIdxs)];
        end_seg = [timestamps(segmentEndIdxs)];
        
        for s = 1:length(start_seg)
            plot([start_seg(s),start_seg(s)],[-5,6],'m');
        end
        for e = 1:length(end_seg)
            plot([end_seg(e),end_seg(e)],[-5,6],'m');
        end
    else
        figure(10); set(gcf, 'name', 'No-PD Peak-Detection')
        subplot(2, 3, subject-10);
        plot(matrix(:,1),energy);
        hold on;
        %plot(matrix(:,1),matrix(:,3), matrix(:,1),matrix(:,4));
        timestamps = matrix(:,1);
        peak_timestamps = timestamps(peakLocInds);
        plot(peak_timestamps, peak, 'r.');
        title(strcat('kav',id));
        xlabel('time (ms)'); ylabel('acceleration (m/s^2)');
        
        start_seg = [timestamps(segmentStartIdxs)];
        end_seg = [timestamps(segmentEndIdxs)];
        
        for s = 1:length(start_seg)
            plot([start_seg(s),start_seg(s)],[-5,6],'m');
        end
        for e = 1:length(end_seg)
            plot([end_seg(e),end_seg(e)],[-5,6],'m');
        end
        
    end
end
peaknumbers

%%
% plot(data_acc_nts(:, 1), data_acc_nts_b(:, 2)); plot(data_acc_nts(:, 1), data_acc_nts_b(:, 3));
% title("Filtered Accelerometer Data w/ Wavelet");
% xlabel("m/s^2");
% ylabel("Time Stamp (ms)");
% legend('x', 'y', 'z');
% hold off;

%% Segmentation
[peaks, peakLocs] = findpeaks(energy,'minPeakHeight',0.02,'minPeakDistance',80);
segmentA = 60;
segmentB = 80;
segmentStartIdxs = peakLocs - segmentA;
segmentEndIdxs = peakLocs + segmentB;
%segmentStartings you have the segment indices.

%     %% wavelet transformation - ignore for now 
% %
%     wt = modwt(matrix(:, 2));
%     figure;
%     %analyze individual wavelets from the decomposition
%     for a  = 1:length(wt(:, 1)) 
%         
% %         subplot(6, 2, a); %m x n plot
% %         plot(wt(a, :));
%     end
%     % take out columns 7-11 for reconstruction
%     wtrec = zeros(size(wt));
%     wtrec(3:6, :) = wt(3:6, :);
%     modified_signal = imodwt(wtrec);    
%     if id(4) == 'A'
%         figure(5); subplot(2, 5, subject);
%         plot(matrix(:, 1),modified_signal);
%     else
%         figure(6); subplot(2, 3, subject-10);
%         plot(matrix(:, 1),modified_signal);
%     end


