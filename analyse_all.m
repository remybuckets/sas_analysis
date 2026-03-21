% analyse_all.m
% TEAM AO
%
%   Navid Khojasteh-Abbasi
%   Jeremy Holguin
%   Zoe Wills
%   Anne Mu
%   Sally Oppenheim
%   Sufiya Al Assad   
%
% Fully automated analysis pipeline for vowel phonation data.
% Processes multiple subjects, both vowels, both microphone channels.
%
% Analyses performed:
%   a. Time-domain analysis (waveforms, extracted windows, detector internals)
%   b. Signal size and noise (RMS, energy, SNR proxy) - both channels
%   c. Frequency-domain analysis (magnitude spectra, essential bandwidth)
%   d. Downsampling and filtering investigation
%   e. Comparisons: vowel vs vowel, channel vs channel, rep vs rep, subject vs subject
%
% All figures saved to ./figures/ as 300 dpi PNGs.
% All metrics saved to all_metrics.mat.
% Summary table printed to console and saved to summary_table.csv.
%
% Requires: Signal Processing Toolbox (butter, filtfilt, freqz)

clear; clc; close all;

% ── Configuration ─────────────────────────────────────────────────────────
SUBJECT_FILES = {'subject01.mat', 'subject02.mat', 'subject03.mat'};
SUBJECT_LABELS = {'S1', 'S2', 'S3'};
FIG_DIR       = 'figures';
BW_FRACTION   = 0.95;    % fraction of energy for essential bandwidth
DS_FACTORS    = [2, 4];  % downsampling factors: Fs/2 and Fs/4
FILT_ORDER    = 6;       % Butterworth filter order
THR_FRAC      = 0.02;    % detector default threshold fraction (for captions only)

if ~exist(FIG_DIR, 'dir'), mkdir(FIG_DIR); end

vowels     = {'aaa', 'eee'};
nSubjects  = numel(SUBJECT_FILES);
nVowels    = 2;
nChannels  = 2;

% all_metrics(s,v) is a struct storing per-rep data for subject s, vowel v
all_metrics = cell(nSubjects, nVowels);

% ══════════════════════════════════════════════════════════════════════════
% MAIN LOOP - iterate over subjects
% ══════════════════════════════════════════════════════════════════════════

for s = 1:nSubjects
    subFile  = SUBJECT_FILES{s};
    subLabel = SUBJECT_LABELS{s};

    if ~exist(subFile, 'file')
        warning('File not found: %s - skipping.', subFile);
        continue
    end

    fprintf('\n%s\n', repmat('=',1,70));
    fprintf('SUBJECT: %s (%s)\n', subLabel, subFile);
    fprintf('%s\n', repmat('=',1,70));

    data    = load(subFile);
    Fs      = double(data.Fs);
    signals = {data.aaa, data.eee};

    % ── Detect repetitions for each vowel ─────────────────────────────────
    opts     = struct();
    D_all    = cell(1, nVowels);
    wins_all = cell(1, nVowels);

    for v = 1:nVowels
        vname = vowels{v};
        X     = signals{v};
        D     = repDetector_x4_strongChannel(X, Fs, opts);
        nReps = numel(D.starts);

        fprintf('\n  [%s | %s] %d rep(s) detected | ch%d used | thr=%.4e\n', ...
                subLabel, vname, nReps, D.chUsed, D.thr);

        if nReps == 0
            warning('  No bursts detected for %s %s - skipping.', subLabel, vname);
            D_all{v}    = D;
            wins_all{v} = {};
            continue
        end
        if nReps < 4
            fprintf('  NOTE: Only %d/4 reps available for %s %s.\n', nReps, subLabel, vname);
        end

        wins = cell(nReps, 1);
        for r = 1:nReps
            wins{r} = X(D.winStarts(r):D.winEnds(r), :);  % keep BOTH channels
        end
        D_all{v}    = D;
        wins_all{v} = wins;
    end

    % ══════════════════════════════════════════════════════════════════════
    % SECTION A - Time-domain analysis
    % ══════════════════════════════════════════════════════════════════════

    for v = 1:nVowels
        vname = vowels{v};
        X     = signals{v};
        D     = D_all{v};
        wins  = wins_all{v};
        nReps = numel(D.starts);
        if nReps == 0, continue; end

        N       = size(X, 1);
        t       = (0:N-1) / Fs;
        colours = lines(nReps);

        % Figure A1 - Both channels full waveform with windows marked
        fig = figure('Visible','off','Position',[100 100 1100 500]);
        for ch = 1:2
            subplot(2,1,ch);
            plot(t, X(:,ch), 'Color', [0.25 0.45 0.75], 'LineWidth', 0.5);
            hold on;
            yl = ylim;
            for r = 1:nReps
                ws = D.winStarts(r)/Fs;
                we = D.winEnds(r)/Fs;
                patch([ws we we ws],[yl(1) yl(1) yl(2) yl(2)], colours(r,:), ...
                      'FaceAlpha',0.15,'EdgeColor',colours(r,:),'LineWidth',1.0);
                if ch == 1
                    text((ws+we)/2, yl(2)*0.82, sprintf('R%d',r), ...
                         'HorizontalAlignment','center','FontSize',8,'Color',colours(r,:));
                end
            end
            if ch == D.chUsed
                ylabel(sprintf('Ch%d amplitude (V)\n[detector channel]', ch));
            else
                ylabel(sprintf('Ch%d amplitude (V)', ch));
            end
            if ch == 2, xlabel('Time (s)'); end
            title(sprintf('%s | "%s" | Channel %d', subLabel, upper(vname), ch));
            grid on;
        end
        sgtitle(sprintf('Figure A1 (%s, %s): Full waveforms - both channels', subLabel, vname));
        caption = sprintf(['Both microphone channels for subject %s, vowel "%s". ' ...
                           'Shaded regions show the %d extracted 0.5 s analysis windows ' ...
                           'centred on each detected phonation burst. ' ...
                           'Channel %d (labelled) was selected by the detector as the stronger channel.'], ...
                           subLabel, vname, nReps, D.chUsed);
        annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('A1_%s_%s_waveforms.png', subLabel, vname)), 'Resolution',300);
        close(fig);

        % Figure A2 - Extracted 0.5 s windows, both channels, all reps
        fig = figure('Visible','off','Position',[100 100 1100 220*nReps]);
        for r = 1:nReps
            for ch = 1:2
                w   = wins{r}(:, ch);
                t_w = (0:numel(w)-1) / Fs;
                subplot(nReps, 2, (r-1)*2 + ch);
                plot(t_w, w, 'Color', colours(r,:), 'LineWidth', 0.8);
                ylabel('Amplitude (V)');
                title(sprintf('Rep %d | Ch%d  [%.3f-%.3f s]', r, ch, D.winStarts(r)/Fs, D.winEnds(r)/Fs));
                grid on;
                if r == nReps, xlabel('Time within window (s)'); end
            end
        end
        sgtitle(sprintf('Figure A2 (%s, %s): Extracted 0.5 s windows - both channels', subLabel, vname));
        caption = sprintf(['The %d extracted 0.5 s analysis windows for subject %s, vowel "%s", ' ...
                           'shown for both microphone channels. Each window is centred on the ' ...
                           'midpoint of the detected phonation burst.'], nReps, subLabel, vname);
        annotation('textbox',[0 0 1 0.03],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('A2_%s_%s_windows.png', subLabel, vname)), 'Resolution',300);
        close(fig);

        % Figure A3 - Detector internals
        fig = figure('Visible','off','Position',[100 100 1000 480]);
        subplot(2,1,1);
        plot(t, D.activity, 'k', 'LineWidth', 0.7); hold on;
        yline(D.thr,'r--','LineWidth',1.5,'DisplayName',sprintf('Threshold (%.3e)',D.thr));
        ylabel('Activity (x^4 smoothed)'); xlabel('Time (s)');
        title(sprintf('%s | "%s" - Energy activity signal', subLabel, vname));
        legend; grid on;
        subplot(2,1,2);
        area(t, double(D.mask),'FaceColor',[0.4 0.7 0.4],'EdgeColor','none','FaceAlpha',0.6);
        ylim([-0.05 1.15]);
        ylabel('Mask'); xlabel('Time (s)');
        title('Binary phonation mask after debounce');
        grid on;
        sgtitle(sprintf('Figure A3 (%s, %s): Detector internals', subLabel, vname));
        caption = sprintf(['Top: x^4 smoothed energy activity with threshold (red dashed, ' ...
                           '%.0f%% of peak). Bottom: binary phonation mask after removal ' ...
                           'of bursts <30 ms and filling of gaps <20 ms.'], THR_FRAC*100);
        annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('A3_%s_%s_detector.png', subLabel, vname)), 'Resolution',300);
        close(fig);
    end

    % ══════════════════════════════════════════════════════════════════════
    % SECTION B - Signal size and noise (BOTH channels)
    % ══════════════════════════════════════════════════════════════════════

    fprintf('\n  -- Section B: RMS / Energy / SNR --\n');

    for v = 1:nVowels
        vname = vowels{v};
        wins  = wins_all{v};
        D     = D_all{v};
        nReps = numel(wins);
        if nReps == 0, continue; end

        rms_vals    = zeros(nReps, nChannels);
        energy_vals = zeros(nReps, nChannels);
        snr_vals    = zeros(nReps, nChannels);

        for r = 1:nReps
            for ch = 1:nChannels
                sig = wins{r}(:, ch);

                rms_vals(r,ch)    = rms(sig);
                energy_vals(r,ch) = sum(sig.^2);
        % Use silence OUTSIDE the burst (before/after) as noise reference
        % instead of window edges (which contain vowel signal)
        
        X_full = signals{v};  % Full signal for this vowel
        N_full = size(X_full, 1);
        
        burst_start = D.starts(r);
        burst_end   = D.ends(r);
        
        % Get 100ms of silence before and after the burst
        pre_len  = round(0.1 * Fs);
        post_len = round(0.1 * Fs);
        
        pre_start  = max(1, burst_start - pre_len);
        pre_end    = max(1, burst_start - 1);
        post_start = min(N_full, burst_end + 1);
        post_end   = min(N_full, burst_end + post_len);
        
        % Extract noise samples from silent regions
        if pre_end >= pre_start && post_end >= post_start
            noise_pre  = X_full(pre_start:pre_end, ch);
            noise_post = X_full(post_start:post_end, ch);
            noise = [noise_pre; noise_post];
        elseif pre_end >= pre_start
            noise = X_full(pre_start:pre_end, ch);
        elseif post_end >= post_start
            noise = X_full(post_start:post_end, ch);
        else
            noise = [];
        end
        
        % Compute SNR using central 80% of window as signal
        n     = numel(sig);
        n80s  = round(0.10*n);
        n80e  = round(0.90*n);
        s_mid = sig(n80s:n80e);
        p_sig = mean(s_mid.^2);
        
        if ~isempty(noise) && numel(noise) > 10
            p_noise = mean(noise.^2);
            if p_noise > 0
                snr_vals(r,ch) = 10*log10(p_sig / p_noise);
            else
                snr_vals(r,ch) = Inf;
            end
        else
            snr_vals(r,ch) = NaN;  
        end


                fprintf('  [%s|%s] Rep%d Ch%d | RMS=%.5f | E=%.5f | SNR=%.2f dB\n', ...
                        subLabel, vname, r, ch, rms_vals(r,ch), energy_vals(r,ch), snr_vals(r,ch));
            end
        end

        % Figure B1 - RMS, energy, and SNR per rep, both channels
        fig = figure('Visible','off','Position',[100 100 1200 380]);
        subplot(1,3,1);
        bar(rms_vals);
        legend('Channel 1','Channel 2','Location','northeast');
        xlabel('Repetition'); ylabel('RMS Amplitude (V)');
        title(sprintf('%s | "%s" - RMS per rep', subLabel, vname));
        xticks(1:nReps); grid on;
        subplot(1,3,2);
        bar(energy_vals);
        legend('Channel 1','Channel 2','Location','northeast');
        xlabel('Repetition'); ylabel('Energy (V^2)');
        title(sprintf('%s | "%s" - Energy per rep', subLabel, vname));
        xticks(1:nReps); grid on;
        subplot(1,3,3);
        bar(snr_vals);
        legend('Channel 1','Channel 2','Location','northeast');
        xlabel('Repetition'); ylabel('SNR proxy (dB)');
        title(sprintf('%s | "%s" - SNR per rep', subLabel, vname));
        xticks(1:nReps); grid on;
        sgtitle(sprintf('Figure B1 (%s, %s): Signal size - both channels', subLabel, vname));
        caption = sprintf(['RMS amplitude (left), signal energy (centre), and SNR proxy (right) ' ...
                           'per repetition for both channels, subject %s, vowel "%s". ' ...
                           'Energy is the sum of squared samples. SNR is the ratio of power ' ...
                           'in the central 80%% of the window to power in the edge 5%% ' ...
                           '(pre/post-burst noise floor).'], subLabel, vname);
        annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('B1_%s_%s_rms_energy_snr.png', subLabel, vname)), 'Resolution',300);
        close(fig);

        % Store in all_metrics
        all_metrics{s,v}.rms    = rms_vals;
        all_metrics{s,v}.energy = energy_vals;
        all_metrics{s,v}.snr    = snr_vals;
        all_metrics{s,v}.nReps  = nReps;
        all_metrics{s,v}.chUsed = D.chUsed;
    end

    % ══════════════════════════════════════════════════════════════════════
    % SECTION C - Frequency-domain analysis (BOTH channels)
    % ══════════════════════════════════════════════════════════════════════

    fprintf('\n  -- Section C: Spectra / Bandwidth --\n');

    for v = 1:nVowels
        vname = vowels{v};
        wins  = wins_all{v};
        D     = D_all{v};
        nReps = numel(wins);
        if nReps == 0, continue; end

        bw_vals = zeros(nReps, nChannels);

        % Figure C1 - Magnitude spectra, both channels overlaid
        fig = figure('Visible','off','Position',[100 100 1100 420]);
        ch_colours = {lines(nReps), lines(nReps)};
        for ch = 1:nChannels
            subplot(1,2,ch);
            legend_str = cell(nReps,1);
            for r = 1:nReps
                sig   = wins{r}(:, ch);
                N     = numel(sig);
                win   = hann(N);
                X_f   = fft(sig .* win, N);
                mag   = abs(X_f(1:floor(N/2)+1));
                mag   = mag / max(mag);
                freqs = (0:floor(N/2)) * Fs / N;

                pow     = mag.^2;
                cum_pow = cumsum(pow) / sum(pow);
                bw_idx  = find(cum_pow >= BW_FRACTION, 1);
                bw_vals(r,ch) = freqs(bw_idx);

                c = ch_colours{ch}(r,:);
                semilogy(freqs/1000, mag, 'Color', c, 'LineWidth', 1.0);
                hold on;
                legend_str{r} = sprintf('Rep%d BW=%.0fHz', r, bw_vals(r,ch));
            end
            xlim([0 Fs/2/1000]);
            ylabel('Normalised Magnitude'); 
            title(sprintf('Ch%d', ch));
            legend(legend_str,'Location','northeast','FontSize',7);
            grid on;
            xlabel('Frequency (kHz)');
        end
        sgtitle(sprintf('Figure C1 (%s, %s): Magnitude spectra - both channels', subLabel, vname));
        caption = sprintf(['Normalised Hann-windowed FFT magnitude spectra for subject %s, ' ...
                           'vowel "%s", all repetitions, both channels. Legend shows the %d%% ' ...
                           'essential bandwidth per repetition.'], subLabel, vname, round(BW_FRACTION*100));
        annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('C1_%s_%s_spectra.png', subLabel, vname)), 'Resolution',300);
        close(fig);

        all_metrics{s,v}.bandwidth = bw_vals;
        fprintf('  [%s|%s] BW Ch1: %s Hz\n', subLabel, vname, num2str(bw_vals(:,1)', '%.0f '));
        fprintf('  [%s|%s] BW Ch2: %s Hz\n', subLabel, vname, num2str(bw_vals(:,2)', '%.0f '));
    end

    % Figure C2 - Spectrograms: aaa vs eee, strong channel
    fig = figure('Visible','off','Position',[100 100 1100 420]);
    for v = 1:nVowels
        vname = vowels{v};
        wins  = wins_all{v};
        D     = D_all{v};
        if isempty(wins), continue; end
        sig = wins{1}(:, D.chUsed);
        subplot(1,2,v);
        spectrogram(sig, hann(256), 200, 512, Fs, 'yaxis');
        ylim([0 8]);
        title(sprintf('%s | "%s" Rep 1 (Ch%d)', subLabel, vname, D.chUsed));
        colorbar;
    end
    sgtitle(sprintf('Figure C2 (%s): Spectrograms - "aaa" vs "eee"', subLabel));
    caption = sprintf(['STFT spectrograms (256-sample Hann window, 200-sample overlap) ' ...
                       'of the first repetition for each vowel, subject %s. ' ...
                       'Colour indicates power spectral density (dB/Hz).'], subLabel);
    annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
    exportgraphics(fig, fullfile(FIG_DIR, sprintf('C2_%s_spectrograms.png', subLabel)), 'Resolution',300);
    close(fig);

    % ══════════════════════════════════════════════════════════════════════
    % SECTION D - Downsampling and filtering
    % ══════════════════════════════════════════════════════════════════════

    fprintf('\n  -- Section D: Downsampling / Filtering --\n');

    for v = 1:nVowels
        vname = vowels{v};
        wins  = wins_all{v};
        D     = D_all{v};
        if isempty(wins), continue; end

        sig_orig  = wins{1}(:, D.chUsed);
        bw_median = median(all_metrics{s,v}.bandwidth(:, D.chUsed));

        % Figure D1 - Spectra at original and downsampled rates
        fig = figure('Visible','off','Position',[100 100 1100 380]);
        subplot(1, numel(DS_FACTORS)+1, 1);
        N     = numel(sig_orig);
        X_f   = fft(sig_orig .* hann(N), N);
        mag   = abs(X_f(1:floor(N/2)+1));
        mag   = mag / max(mag);
        freqs = (0:floor(N/2)) * Fs / N;
        semilogy(freqs/1000, mag, 'Color',[0.2 0.4 0.8], 'LineWidth',0.9);
        xlim([0 Fs/2/1000]); xlabel('Freq (kHz)'); ylabel('Norm. Mag.');
        title(sprintf('Original\nFs=%d Hz', Fs)); grid on;

        for di = 1:numel(DS_FACTORS)
            fac   = DS_FACTORS(di);
            Fs_ds = Fs / fac;
            fc_aa = (Fs_ds/2) * 0.9;
            [b_aa, a_aa] = butter(FILT_ORDER, fc_aa/(Fs/2), 'low');
            sig_ds = filtfilt(b_aa, a_aa, sig_orig);
            sig_ds = sig_ds(1:fac:end);

            subplot(1, numel(DS_FACTORS)+1, di+1);
            N_ds   = numel(sig_ds);
            X_ds   = fft(sig_ds .* hann(N_ds), N_ds);
            mag_ds = abs(X_ds(1:floor(N_ds/2)+1));
            mag_ds = mag_ds / max(mag_ds);
            fr_ds  = (0:floor(N_ds/2)) * Fs_ds / N_ds;
            semilogy(fr_ds/1000, mag_ds, 'Color',[0.8 0.3 0.2], 'LineWidth',0.9);
            xlim([0 Fs_ds/2/1000]); xlabel('Freq (kHz)'); ylabel('Norm. Mag.');
            title(sprintf('Fs/%d\n=%d Hz', fac, Fs_ds)); grid on;

            fprintf('  [%s|%s] DS x%d: Fs=%d Hz | RMS=%.5f\n', ...
                    subLabel, vname, fac, Fs_ds, rms(sig_ds));
        end
        sgtitle(sprintf('Figure D1 (%s, %s): Effect of downsampling on spectrum', subLabel, vname));
        caption = sprintf(['Normalised spectra at original and downsampled rates for %s "%s" Rep 1. ' ...
                           'A %d-th order Butterworth anti-aliasing LPF (fc=90%% new Nyquist) ' ...
                           'is applied before each downsampling step.'], subLabel, vname, FILT_ORDER);
        annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('D1_%s_%s_downsampling.png', subLabel, vname)), 'Resolution',300);
        close(fig);

        % Figure D2 - Designed LPF response and filtered waveform
        fc_design = min(bw_median * 1.10, Fs/2 * 0.90);
        [b_lp, a_lp] = butter(FILT_ORDER, fc_design/(Fs/2), 'low');
        sig_filtered  = filtfilt(b_lp, a_lp, sig_orig);

        fig = figure('Visible','off','Position',[100 100 950 400]);
        subplot(1,2,1);
        freqz(b_lp, a_lp, 1024, Fs);
        title(sprintf('%s | "%s" LPF (fc=%.0fHz)', subLabel, vname, fc_design));
        subplot(1,2,2);
        t_w = (0:numel(sig_orig)-1)/Fs;
        plot(t_w, sig_orig,    'Color',[0.7 0.7 0.7],'LineWidth',0.8,'DisplayName','Original');
        hold on;
        plot(t_w, sig_filtered,'Color',[0.2 0.6 0.3],'LineWidth',1.0,'DisplayName','Filtered');
        xlabel('Time (s)'); ylabel('Amplitude (V)');
        title('Original vs filtered (Rep 1)'); legend; grid on;
        sgtitle(sprintf('Figure D2 (%s, %s): Designed lowpass filter', subLabel, vname));
        caption = sprintf(['Left: %d-th order Butterworth LPF (fc=%.0f Hz, 10%% above ' ...
                           'median %d%% BW of %.0f Hz). Right: filtered vs original ' ...
                           'waveform for %s "%s" Rep 1.'], ...
                           FILT_ORDER, fc_design, round(BW_FRACTION*100), bw_median, subLabel, vname);
        annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('D2_%s_%s_filter.png', subLabel, vname)), 'Resolution',300);
        close(fig);

        % Figure D3 - Metric change: original, designed LPF, Fs/2, Fs/4
        all_configs = {'Original', sprintf('LPF\n%.0fHz',fc_design), 'Fs/2', 'Fs/4'};
        rms_cmp = zeros(1,4);
        bw_cmp  = zeros(1,4);

        % baseline: original signal
        rms_cmp(1) = rms(sig_orig);
        N_orig = numel(sig_orig);
        X_orig = fft(sig_orig .* hann(N_orig), N_orig);
        mag_orig = abs(X_orig(1:floor(N_orig/2)+1));
        fr_orig  = (0:floor(N_orig/2)) * Fs / N_orig;
        cum_orig = cumsum(mag_orig.^2) / sum(mag_orig.^2);
        bw_cmp(1) = fr_orig(find(cum_orig >= BW_FRACTION, 1));

        % designed LPF at original Fs (no downsampling)
        rms_cmp(2) = rms(sig_filtered);
        X_lp   = fft(sig_filtered .* hann(N_orig), N_orig);
        mag_lp = abs(X_lp(1:floor(N_orig/2)+1));
        cum_lp = cumsum(mag_lp.^2) / sum(mag_lp.^2);
        bw_cmp(2) = fr_orig(find(cum_lp >= BW_FRACTION, 1));

        % downsampled versions
        for di = 1:numel(DS_FACTORS)
            fac   = DS_FACTORS(di);
            Fs_ds = Fs / fac;
            fc_aa = (Fs_ds/2) * 0.9;
            [b_aa, a_aa] = butter(FILT_ORDER, fc_aa/(Fs/2), 'low');
            sig_ds = filtfilt(b_aa, a_aa, sig_orig);
            sig_ds = sig_ds(1:fac:end);
            rms_cmp(di+2) = rms(sig_ds);
            N_ds   = numel(sig_ds);
            X_ds   = fft(sig_ds .* hann(N_ds), N_ds);
            mag_ds = abs(X_ds(1:floor(N_ds/2)+1));
            fr_ds  = (0:floor(N_ds/2)) * Fs_ds / N_ds;
            cum_ds = cumsum(mag_ds.^2) / sum(mag_ds.^2);
            bw_cmp(di+2) = fr_ds(find(cum_ds >= BW_FRACTION, 1));
        end

        fig = figure('Visible','off','Position',[100 100 900 350]);
        bar_colours = [0.3 0.6 0.9; 0.2 0.7 0.4; 0.8 0.5 0.2; 0.75 0.3 0.3];
        subplot(1,2,1);
        b = bar(rms_cmp, 'FaceColor','flat');
        for ci = 1:4, b.CData(ci,:) = bar_colours(ci,:); end
        set(gca,'XTickLabel',all_configs);
        ylabel('RMS (V)');
        title(sprintf('%s | "%s" - RMS', subLabel, vname)); grid on;
        subplot(1,2,2);
        b = bar(bw_cmp, 'FaceColor','flat');
        for ci = 1:4, b.CData(ci,:) = bar_colours(ci,:); end
        set(gca,'XTickLabel',all_configs);
        ylabel(sprintf('BW_{%d%%} (Hz)', round(BW_FRACTION*100)));
        title(sprintf('%s | "%s" - BW', subLabel, vname)); grid on;
        sgtitle(sprintf('Figure D3 (%s, %s): Effect of filtering and downsampling on metrics', subLabel, vname));
        caption = sprintf(['RMS (left) and %d%% essential bandwidth (right) for %s "%s" Rep 1 ' ...
                           'under four conditions: original signal, designed LPF at full Fs ' ...
                           '(fc=%.0f Hz), and downsampled to Fs/2 and Fs/4 with anti-aliasing.'], ...
                           round(BW_FRACTION*100), subLabel, vname, fc_design);
        annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
                    'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
        exportgraphics(fig, fullfile(FIG_DIR, sprintf('D3_%s_%s_metrics_vs_fs.png', subLabel, vname)), 'Resolution',300);
        close(fig);
    end

end % end subject loop

% ══════════════════════════════════════════════════════════════════════════
% SECTION E - Cross-subject and cross-condition comparisons
% ══════════════════════════════════════════════════════════════════════════

fprintf('\n%s\n', repmat('=',1,70));
fprintf('SECTION E - Comparisons\n');
fprintf('%s\n', repmat('=',1,70));

% ── E1: Vowel comparison (aaa vs eee) - RMS and BW, strong channel ────────
% Collect mean RMS and BW per subject per vowel
mean_rms = nan(nSubjects, nVowels, nChannels);
mean_bw  = nan(nSubjects, nVowels, nChannels);
mean_snr = nan(nSubjects, nVowels, nChannels);

for s = 1:nSubjects
    for v = 1:nVowels
        m = all_metrics{s,v};
        if isempty(m), continue; end
        for ch = 1:nChannels
            mean_rms(s,v,ch) = mean(m.rms(:,ch));
            mean_bw(s,v,ch)  = mean(m.bandwidth(:,ch));
            mean_snr(s,v,ch) = mean(m.snr(:,ch));
        end
    end
end

% Figure E1 - Vowel comparison: aaa vs eee (grouped by subject)
fig = figure('Visible','off','Position',[100 100 1100 380]);
metric_names = {'Mean RMS (V)', sprintf('Mean BW_{%d%%} (Hz)', round(BW_FRACTION*100)), 'Mean SNR (dB)'};
metric_data  = {mean_rms, mean_bw, mean_snr};
ch_plot      = 1;   % plot strong channel (ch1); ch2 shown in E2

for mi = 1:3
    subplot(1,3,mi);
    dat = squeeze(metric_data{mi}(:,:,ch_plot));   % [nSubjects x nVowels]
    bar(dat);
    set(gca,'XTickLabel', SUBJECT_LABELS);
    legend('aaa','eee','Location','northeast');
    ylabel(metric_names{mi});
    xlabel('Subject');
    title(metric_names{mi});
    grid on;
end
sgtitle('Figure E1: Vowel comparison (aaa vs eee) - Channel 1');
caption = ['Mean RMS, essential bandwidth, and SNR proxy for vowels "aaa" and "eee" ' ...
           'across all subjects (Channel 1). Each group of bars represents one subject.'];
annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
            'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
exportgraphics(fig, fullfile(FIG_DIR, 'E1_vowel_comparison.png'), 'Resolution',300);
close(fig);

% Figure E2 - Channel comparison (Ch1 vs Ch2) per vowel
fig = figure('Visible','off','Position',[100 100 1100 420]);
for v = 1:nVowels
    vname = vowels{v};
    subplot(1,2,v);
    dat_ch1 = squeeze(mean_rms(:,v,1));
    dat_ch2 = squeeze(mean_rms(:,v,2));
    bar([dat_ch1, dat_ch2]);
    set(gca,'XTickLabel', SUBJECT_LABELS);
    legend('Channel 1','Channel 2','Location','northeast');
    ylabel('Mean RMS (V)');
    xlabel('Subject');
    title(sprintf('"%s" - Ch1 vs Ch2 RMS', vname));
    grid on;
end
sgtitle('Figure E2: Channel comparison - mean RMS per subject per vowel');
caption = ['Mean RMS amplitude for Channel 1 and Channel 2 across subjects and vowels. ' ...
           'Differences between channels reflect microphone placement and sensitivity.'];
annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
            'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
exportgraphics(fig, fullfile(FIG_DIR, 'E2_channel_comparison.png'), 'Resolution',300);
close(fig);

% Figure E3 - Subject comparison: all subjects, both vowels, RMS
fig = figure('Visible','off','Position',[100 100 1000 400]);
subplot(1,2,1);
bar(squeeze(mean_rms(:,:,1)));
set(gca,'XTickLabel', SUBJECT_LABELS);
legend('aaa','eee','Location','northeast');
ylabel('Mean RMS (V)'); xlabel('Subject');
title('Ch1 - RMS across subjects'); grid on;
subplot(1,2,2);
bar(squeeze(mean_bw(:,:,1)));
set(gca,'XTickLabel', SUBJECT_LABELS);
legend('aaa','eee','Location','northeast');
ylabel(sprintf('Mean BW_{%d%%} (Hz)', round(BW_FRACTION*100)));
xlabel('Subject');
title(sprintf('Ch1 - BW_{%d%%} across subjects', round(BW_FRACTION*100))); grid on;
sgtitle('Figure E3: Subject comparison - RMS and bandwidth (Channel 1)');
caption = ['Mean RMS amplitude (left) and essential bandwidth (right) for all subjects, ' ...
           'both vowels, Channel 1. Inter-subject variability reflects differences in ' ...
           'vocal characteristics and recording conditions.'];
annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
            'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
exportgraphics(fig, fullfile(FIG_DIR, 'E3_subject_comparison.png'), 'Resolution',300);
close(fig);

% Figure E4 - Repetition consistency: RMS per rep across subjects, per vowel
fig = figure('Visible','off','Position',[100 100 1100 420]);
for v = 1:nVowels
    vname = vowels{v};
    subplot(1,2,v);
    hold on;
    subj_colours = lines(nSubjects);
    for s = 1:nSubjects
        m = all_metrics{s,v};
        if isempty(m), continue; end
        nR = m.nReps;
        plot(1:nR, m.rms(:,1), '-o', 'Color', subj_colours(s,:), ...
             'LineWidth',1.2,'DisplayName', SUBJECT_LABELS{s});
    end
    xlabel('Repetition'); ylabel('RMS (V)');
    title(sprintf('"%s" - RMS per rep (Ch1)', vname));
    legend('Location','best','FontSize',8); grid on;
    xlim([0.5 4.5]); xticks(1:4);
end
sgtitle('Figure E4: Repetition consistency - RMS per repetition (Channel 1)');
caption = ['RMS amplitude per repetition for each subject (Channel 1). ' ...
           'Stability across repetitions indicates consistent phonation; ' ...
           'variation may reflect fatigue, pitch drift, or recording artefacts.'];
annotation('textbox',[0 0 1 0.04],'String',caption,'EdgeColor','none', ...
            'FontSize',7,'FontAngle','italic','HorizontalAlignment','center');
exportgraphics(fig, fullfile(FIG_DIR, 'E4_repetition_consistency.png'), 'Resolution',300);
close(fig);

% ── Summary table ──────────────────────────────────────────────────────────
fprintf('\n%-4s %-5s %-4s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
        'Subj','Vowel','Ch','RMS mean','RMS std','BW mean','BW std','SNR mean','SNR std');
fprintf('%s\n', repmat('-',1,80));

fid = fopen('summary_table.csv','w');
fprintf(fid,'Subject,Vowel,Channel,RMS_mean,RMS_std,BW_mean,BW_std,SNR_mean,SNR_std\n');

for s = 1:nSubjects
    for v = 1:nVowels
        vname = vowels{v};
        m = all_metrics{s,v};
        if isempty(m), continue; end
        for ch = 1:nChannels
            rms_m = mean(m.rms(:,ch));
            rms_s = std(m.rms(:,ch));
            bw_m  = mean(m.bandwidth(:,ch));
            bw_s  = std(m.bandwidth(:,ch));
            snr_m = mean(m.snr(:,ch));
            snr_s = std(m.snr(:,ch));

            fprintf('%-4s %-5s %-4d %-10.5f %-10.5f %-10.1f %-10.1f %-10.2f %-10.2f\n', ...
                    SUBJECT_LABELS{s}, vname, ch, rms_m, rms_s, bw_m, bw_s, snr_m, snr_s);
            fprintf(fid,'%s,%s,%d,%.5f,%.5f,%.1f,%.1f,%.2f,%.2f\n', ...
                    SUBJECT_LABELS{s}, vname, ch, rms_m, rms_s, bw_m, bw_s, snr_m, snr_s);
        end
    end
end
fclose(fid);

% ── Save all metrics ────────────────────────────────────────────────────────
save('all_metrics.mat', 'all_metrics', 'vowels', 'SUBJECT_LABELS', 'Fs', 'BW_FRACTION');
fprintf('\nAll metrics saved to all_metrics.mat\n');
fprintf('Summary table saved to summary_table.csv\n');
fprintf('All figures saved to ./%s/\n', FIG_DIR);
