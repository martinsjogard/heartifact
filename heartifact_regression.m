addpath('/cluster/manoach/toolboxes/homemade_matlab_functions')

dataname = 'WFX002_DREEM_1';
winsize = 1000;
filename = ['heartifact_detection_' dataname ...
            '_win' num2str(winsize) 's_'];
load(filename);

edf_file = ['/autofs/cluster/manoach/MartinS/Seagate_iEEG/FraCluster_copy/' dataname '.edf'];
[hdr, data] = EDFread(edf_file);

EEGchs = {'EEGF7O1','EEGF7O2','EEGF8O1','EEGF8O2'};
[EEGchs, ia, ~] = intersect(EEGchs, hdr.label);
data = data(ia, :);
hdr.label = hdr.label(ia);

Fs = hdr.frequency(1);
datach = data(1,:);

template = nanmean(handles.savewaves,1);
Fs = handles.fs;
locs = handles.locs;

eeg_raw = datach;

% Build artifact timecourse (same length as eeg_raw)
artifact_model = zeros(size(eeg_raw));
for i = 1:length(locs)
    idx_start = locs(i) - floor(size(handles.savewaves,2)/2);
    idx_end = locs(i) + floor(size(handles.savewaves,2)/2);
    if idx_start > 0 && idx_end <= length(eeg_raw)
        artifact_model(idx_start:idx_end) = artifact_model(idx_start:idx_end) + template;
    end
end

t = (1:numel(eeg_raw))/Fs;
usetime = 1:numel(eeg_raw);


% Compute scaling factor and regress out artifact
X = artifact_model(:);
Y = eeg_raw(:);
beta = (X' * X) \ (X' * Y);  % least-squares solution
artifact_est = beta * artifact_model;
eeg_cleaned = eeg_raw - artifact_est;

[b, a] = butter(4, 5 / (Fs/2), 'high');  % 4th-order IIR highpass at 0.5 Hz
eeg_filtered = filtfilt(b, a, eeg_raw);
eeg_cleaned_filtered = filtfilt(b, a, eeg_segment);

for ii = 1:10
    R=randperm(numel(eeg_raw)-10000);
    R=R(1);
    plot_xlim = [t(R) t(R+999)];   

    figure;hold on;

    subplot(3,3,1);
    plot(template)
    hold on;
    %plot(beta*template,'--')
    ylim([min(template)-nanstd(template) ...
        max(template) + nanstd(template)])
    title('Template estimate')

    % cleaned_template = nan(size(handles.savewaves));
    % l1 = locs-(size(handles.savewaves,2)-1)/2;
    % l2 = locs+(size(handles.savewaves,2)-1)/2;
    % for iih = 1:numel(L_ind)     
    %     cleaned_template(iih,:) = eeg_cleaned(l1(iih):l2(iih));
    % end
    % cleaned_template = nanmean(cleaned_template,1);

    L_ind = find((locs>R) & (locs<R+999));
    localtemplate = nanmean(handles.savewaves(L_ind,:),1);
    if sum(~isnan(localtemplate))>numel(localtemplate)/4
        subplot(3,3,4);
        plot(localtemplate)
        hold on;
        ylim([min(localtemplate)-nanstd(localtemplate) ...
            max(localtemplate) + nanstd(localtemplate)])
        title('Local detection average');
    end


    subplot(3,3,[2:3])
    plot(t,eeg_raw);
    xticklabels(t)
    hold on;
    plot(t(locs), eeg_raw(locs), 'ro');
    xlim(plot_xlim);
    title(['Raw EEG (' num2str(numel(L_ind)*6) ' bpm)'])
    
    
    subplot(3,3,[5:6]);
    % plot(t,eeg_filtered);
    % xticklabels(t)
    % hold on;
    % plot(t(locs), eeg_filtered(locs), 'ro');
    % xlim(plot_xlim);
    % title(['Highpass EEG'])
    plot(t,beta*artifact_model);
    xticklabels(t)
    xlim(plot_xlim);
    title('Artifact model')
    
    
    subplot(3,3,[8:9]);
    plot(t,eeg_cleaned);
    xticklabels(t)
    hold on;
    plot(t(locs), eeg_raw(locs), 'ro');
    hold on;
    plot(t(locs), eeg_cleaned(locs), 'gx');
    xlim(plot_xlim);
    title('Cleaned EEG')
    
    sgtitle('Heartbeat Artifact Removal via Template Regression');% Set paper size and figure dimensions
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 12 8]);     
    set(gcf, 'PaperSize', [12 8]); 
    
    
    saveas(gcf,[filename num2str(ii) '.pdf'])

    close all
end

