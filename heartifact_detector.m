clear, clc, close all;
addpath('/cluster/manoach/toolboxes/homemade_matlab_functions')

do_plot_figures = false;

% Load data
subname = 'WAS056_DREEM_2';

for chii = 1:4
    edf_file = ['/autofs/cluster/manoach/MartinS/Seagate_iEEG/FraCluster_copy/' subname '.edf'];
    use_ch = chii;
    [hdr, data] = EDFread(edf_file);
    
    EEGchs = {'EEGF7O1','EEGF7O2','EEGF8O1','EEGF8O2'};
    [EEGchs, ia, ~] = intersect(EEGchs, hdr.label);
    data = data(ia, :);
    hdr.label = hdr.label(ia);
    
    Fs = hdr.frequency(1);
    datach = data(use_ch,:);
    
    % Parameters
    segmsize = 1000;
    
    handles.fs = Fs;
    handles.datach = datach;
    handles.window_size = segmsize * Fs;
    handles.step_size = handles.window_size;
    handles.start_sample = 1;
    handles.tfr_wind_length = round(handles.window_size / 20);
    handles.tfr_wind_overlap = round(handles.tfr_wind_length * 0.8);
    handles.subname = subname;
    handles.savewaves = [];
    
    
    nWins = ceil(numel(datach)/handles.window_size);
    
    for winii = 1:nWins
        disp(['** Window ' num2str(winii) ' / ' num2str(nWins) ' **'])
        close all
    
        
        if winii ~= 1
            handles.start_sample = handles.start_sample + handles.step_size;
            handles.start_sample = max(1, min(handles.start_sample, ...
                    length(handles.datach) - handles.window_size + 1));
        end
    
        start_sample = handles.start_sample;
        end_sample = min(start_sample + handles.window_size - 1, length(handles.datach));
        eeg_segment = handles.datach(start_sample:end_sample);
        time_axis = (start_sample:end_sample) / handles.fs;
        
        % Compute and plot the power spectrum
        [SP,F] = pwelch(handles.datach',handles.fs,[],[],handles.fs);
    
        eeg_segment_clean = eeg_segment;
        eeg_segment_clean(eeg_segment_clean>2 * mad(eeg_segment))=NaN;
        eeg_segment_clean(eeg_segment_clean<3 * mad(eeg_segment))=NaN;
    
            
        % Detect large peaks ~QRS-like transients
        min_peak_height1 = 1 * std(eeg_segment); % Adjust multiplier if needed
        min_peak_height2 = 1 * mad(eeg_segment); % Adjust multiplier if needed
        
        min_peak_distance = round(0.6 * handles.fs);    % ~60 bpm max heart rate
        
        [peaks1A, locs1A] = findpeaks(-eeg_segment, ...
            'MinPeakHeight', min_peak_height1, ...
            'MinPeakDistance', min_peak_distance);
        
        [peaks2A, locs2A] = findpeaks(-eeg_segment, ...
            'MinPeakHeight', min_peak_height2, ...
            'MinPeakDistance', min_peak_distance);
        
        
        %% FILTER AND RERUN PEAKFINDER
        [b, a] = butter(4, 5 / (handles.fs/2), 'high');  % 4th-order IIR highpass at 0.5 Hz
        eeg_filtered = filtfilt(b, a, eeg_segment);
        
        % figure;hold on;
        % subplot(3,3,[2 3]);
        % plot(time_axis, eeg_segment);
        % 
        % subplot(3,3,[5 6]);
        % plot(time_axis, eeg_filtered);
        
        eeg_filtered_clean = eeg_filtered;
        eeg_filtered_clean(eeg_filtered_clean<nanmedian(eeg_filtered)-(7*mad(eeg_filtered)))=NaN;
        eeg_filtered_clean(eeg_filtered_clean>nanmedian(eeg_filtered)+(2*mad(eeg_filtered)))=NaN;
        
      
        % Detect large peaks ~QRS-like transients
        
        min_peak_height1 = 1 * std(eeg_filtered_clean); % Adjust multiplier if needed
        min_peak_height2 = 1 * mad(eeg_filtered_clean); % Adjust multiplier if needed
        
        min_peak_distance = round(0.6 * handles.fs);    % ~60 bpm max heart rate
        
        [peaks1, locs1] = findpeaks(-eeg_filtered_clean, ...
            'MinPeakHeight', min_peak_height1, ...
            'MinPeakDistance', min_peak_distance);
        
        [peaks2, locs2] = findpeaks(-eeg_filtered_clean, ...
            'MinPeakHeight', min_peak_height2, ...
            'MinPeakDistance', min_peak_distance);
        
        locsX = locs2 + handles.start_sample -1;  
        
        
        %% AVERAGE HEARTBEAT WAVEFORM
            
        locs = locs2;
        eeg_seg = eeg_filtered_clean;
        
        pre_samples = round(0.2 * handles.fs);  % 100 ms before
        post_samples = round(0.2 * handles.fs); % 100 ms after
        epoch_len = pre_samples + post_samples + 1;
        
        % Preallocate
        valid_epochs = [];
        for i = 1:length(locs)
            idx_start = locs(i) - pre_samples;
            idx_end = locs(i) + post_samples;
            if idx_start > 0 && idx_end <= length(eeg_seg)
                segment = eeg_seg(idx_start:idx_end);
                valid_epochs = [valid_epochs; segment];  % each row is one heartbeat epoch
            end
        end
        
        % Compute average
        mean_waveform = nanmean(valid_epochs, 1);
        t = (-pre_samples:post_samples) / handles.fs;  % time axis in seconds
        
        
        
        if do_plot_figures
            figure;
            subplot(3, 7, [1]);
            plot(F,pow2db(SP),'LineWidth',2)
            hold on;
            plot(F(F<3),pow2db(SP(F<3)),'r','LineWidth',2)
            set(gca,'XLim',[0,30])
            
            subplot(3, 7, [8]);
            plot(F,pow2db(SP),'r','LineWidth',2)
            set(gca,'XLim',[0,3])
            subplot(3,7,[2:4]);
            plot(eeg_segment)
            hold on;plot(eeg_segment_clean,'r')
            title('Raw signal');
    
            % Plot STD detections on raw eeg
            subplot(3,7,[9:11]);
            plot(time_axis, eeg_segment); hold on;
            plot(time_axis(locs1A), eeg_segment(locs1A), 'ro');
            title(['Detected Peaks STD: ' num2str(numel(locs1A)*(60/(handles.window_size/handles.fs))) ' bpm'  ]);
            
            % Plot MAD detections on raw eeg
            subplot(3,7,[16:18]);
            plot(time_axis, eeg_segment); hold on;
            plot(time_axis(locs2A), eeg_segment(locs2A), 'ro');
            title(['Detected Peaks MAD: ' num2str(numel(locs2A)*(60/(handles.window_size/handles.fs))) ' bpm'  ]);
            xlabel('Time (s)');
            
            subplot(3,7,[5:7]);
            plot(eeg_filtered)
            hold on;plot(eeg_filtered_clean,'r')
            title('5Hz high-pass filtered signal');
           
            % Plot STD detections on clean filtered eeg
            subplot(3,7,[12:14]);
            plot(time_axis, eeg_filtered_clean); hold on;
            plot(time_axis(locs1), eeg_filtered_clean(locs1), 'ro');
            title(['Detected Peaks STD: ' num2str(numel(locs1)*(60/(handles.window_size/handles.fs))) ' bpm'  ]);
            
            % Plot MAD detections on clean filtered eeg
            subplot(3,7,[19:21]);
            plot(time_axis, eeg_filtered_clean); hold on;
            plot(time_axis(locs2), eeg_filtered_clean(locs2), 'ro');
            title(['Detected Peaks MAD: ' num2str(numel(locs2)*(60/(handles.window_size/handles.fs))) ' bpm'  ]);
            xlabel('Time (s)');
        
            sgtitle(sprintf('Start Sample: %d', start_sample));
            % Set paper size and figure dimensions
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperPosition', [0 0 12 8]);     
            set(gcf, 'PaperSize', [12 8]); 
    
            subplot(3,7,15);hold on;
            for iix = 1:size(valid_epochs,1)
                plot(t,valid_epochs(iix,:),'k','LineWidth', 0.5);
            
            end
            
            plot(t, mean_waveform,'b', 'LineWidth', 2);
            xlabel('Time (s)');
            ylabel('Amplitude (\muV)');
            title(sprintf('HB Waveform (N = %d)', size(valid_epochs, 1)));
            grid on;
        
            output_filename = ['heartbeats_detected_' handles.subname ...
                '_ch' num2str(use_ch) '_win' num2str(handles.window_size/handles.fs) 's_' num2str(handles.start_sample)];
    
            disp('SAVING')
        end
        
        
        if isfield(handles, 'savewaves') && ~isempty(handles.savewaves)
            handles.savewaves = [handles.savewaves; valid_epochs];
        else
            handles.savewaves = valid_epochs;
        end
        
        if isfield(handles, 'locs') && ~isempty(handles.locs)
            handles.locs = [handles.locs locsX];
        else
            handles.locs = locsX;
        end 
        
        
        
        
        % 
        % handlesx = guidata(gcf);
    end
    
    
    all_waves = handles.savewaves;
    
    
    output_filename = ['heartifact_detection_' handles.subname ...
                '_ch' num2str(use_ch) '_win' num2str(handles.window_size/handles.fs) 's_'];
    save(output_filename,'handles');

    clearvars -except do_plot_figures subname chii
end


% all_waves = all_waves(sum(isnan(all_waves),2)<5,:);
% figure;hold on
% for iix = 1:size(all_waves,1)
%     plot(all_waves(iix,:),'k','LineWidth', 0.5);
% end

% plot(nanmean(all_waves,1),'b', 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('Amplitude (\muV)');
% title(sprintf('Average Heartbeat Waveform (N = %d)', size(all_waves, 1)));
% grid on;
% output_filename = ['Avg_heartbeats_all_' handlesx.subname ...
%         '_win' num2str(handlesx.window_size/handles.fs)];
% saveas(gcf,[output_filename '.png']);  % 300 dpi PDF

