%% =======================================================================
%  FULL BP ESTIMATION PIPELINE (UPDATED FOR BP_annotate BEHAVIOR)
%  - Beat-to-beat estimation
%  - Uses only last 15 seconds from each 20-sec window
%  - Uses BP_annotate resampled 200 Hz signal
% =======================================================================

clear; clc;

load paired_data.mat   % contains paired{k}.signal + paired{k}.ref

numSamples = length(paired);
numWindows = 3;
t_orig = 0.06;
fs_orig = 1/t_orig;       % original sampling rate (modify if needed)
fs_annot = 200;       % BP_annotate resampled frequency (fixed)

windowLength = 20;    % seconds
N_win = windowLength * fs_orig;

Results = struct();

all_est_SBP = [];
all_ref_SBP = [];
all_est_DBP = [];
all_ref_DBP = [];


%% =======================================================================
%   LOOP THROUGH ALL SAMPLES
% =======================================================================
store_calib = struct();
num_pulses = 0;
num_segments = 0;
for k = 1:numSamples

    fprintf("Processing sample %d/%d...\n", k, numSamples);
    sig = paired{k}.signal.magnitude*-1; % flipped for all samples    
    time_sig = paired{k}.signal.time;
    ref = paired{k}.ref;
    
    fprintf("Obtaining Calibration Parameters ... \n");
    calibration_params = getcalibrationparams(sig, ref, fs_orig);
    store_calib(k).calib = calibration_params;

    %full signal annotate
    full_annot = struct();
    [full_annot.diaTroughIndex, full_annot.sysPeakIndex, full_annot.notchIndex, full_annot.dicroticIndex, full_annot.resampTime, full_annot.resampSignal]   = BP_annotate(sig, fs_orig,0);
    full_trough_id = full_annot.diaTroughIndex;
    full_peak_id = full_annot.sysPeakIndex;

    for w = 1:numWindows
        fprintf("Starting Sample %d, Window %d ... \n", k, w);
        %% ---- Define 20-sec block ----
        t_start = (w-1)*windowLength;
        t_end   = w*windowLength; 
        
        idx = (time_sig >= t_start) & (time_sig < t_end);
        time_win = time_sig(idx);
        sig_win = sig(idx);

        %% ---- Use full 20 seconds for BP_annotate ----
        skipSec = 0; % Adjust to eliminate speicific time window
        valid_idx = time_sig(idx) >= (t_start + skipSec);
        sig_annot = sig_win(valid_idx);
        
        
        if length(sig_annot) < fs_orig  % must be at least 1 second
            fprintf("Window %d of sample %d too short!\n", w, k);
            continue;
        end
        

        %% ---- Run BP_annotate on truncated 15-sec window ----
        try
            fprintf("Compiling into annot \n")
            annot = struct();
            [annot.diaTroughIndex, annot.sysPeakIndex, annot.notchIndex, annot.dicroticIndex, annot.resampTime, annot.resampSignal]   = BP_annotate(sig_annot, fs_orig,0);
        catch
            annot = [];
        end

        if isempty(annot) || isempty(annot.sysPeakIndex)
            % No beats detected
            SBP_vec = [];
            DBP_vec = [];
            SBP_est_win = NaN;
            DBP_est_win = NaN;
            ref_win = extract_ref_window(ref, t_start, t_end);
        else
            %% ---- Extract beat features from resampled 200 Hz signal ----
            [features_raw, peakTimes, stats] = extract_features(annot);

            % Apply beat-level outlier rejection
            idx_valid = beat_filter(features_raw, stats);
            features  = features_raw(idx_valid, :);
            peakTimes = peakTimes(idx_valid);


            %% ---- Estimate beat-to-beat BP ----
            [SBP_vec, DBP_vec, corrected_sig, corrected_timestamps, rej_pulse, num_seg] = getwindowBPestimates(sig_annot,annot,calibration_params, fs_annot, stats);
            
            num_pulses = num_pulses + length(SBP_vec);
            num_segments = num_segments + num_seg;

            %% ---- Window-level estimated BP ----
            SBP_est_win = mean(SBP_vec);
            DBP_est_win = mean(DBP_vec);

            %% ---- Reference BP for this 20-sec window ----
            ref_win = extract_ref_window(ref, t_start, t_end);
        end

        %% ---- Window-level error ----
        SBP_err = SBP_est_win - ref_win.SBP_mean;
        DBP_err = DBP_est_win - ref_win.DBP_mean;
        num_rej = rej_pulse;

        %% ---- Store results ----
        Results(k).window(w).stats = stats;        % FULL FEATURE SUMMARY 
        Results(k).window(w).features = features;  % beat-level features

        Results(k).window(w).SBP_beat = SBP_vec;
        Results(k).window(w).DBP_beat = DBP_vec;

        Results(k).window(w).num_rej = num_rej;

        Results(k).window(w).SBP_est = SBP_est_win;
        Results(k).window(w).DBP_est = DBP_est_win;

        Results(k).window(w).ref_SBP = ref_win.SBP_mean;
        Results(k).window(w).ref_DBP = ref_win.DBP_mean;

        Results(k).window(w).SBP_err = SBP_err;
        Results(k).window(w).DBP_err = DBP_err;

        Results(k).window(w).peakTimes = peakTimes;

        %% ---- For Bland–Altman ----
        all_est_SBP(end+1) = SBP_est_win;
        all_ref_SBP(end+1) = ref_win.SBP_mean;

        all_est_DBP(end+1) = DBP_est_win;
        all_ref_DBP(end+1) = ref_win.DBP_mean;
        fprintf("End of Sample %d, Window %d ... \n", k, w);

        if k == 25 && w == 2
            raw_sig = paired{k}.signal.magnitude;
            raw_time = paired{k}.signal.time;
            figure
            subplot(3,1,1)
            plot(raw_time, raw_sig)
            title('Raw Signal')
            subplot(3,1,2)
            plot(time_win, sig_annot)
            title('20-Sec Window')
            subplot(3,1,3)
            plot(corrected_timestamps, corrected_sig)
            ylim([40 150])
            title('Corrected Signal')

        end

    end
end


%% =======================================================================
%  PART 1 -- Bland–Altman Plots
% =======================================================================
figure;
bland_altman(all_est_SBP, all_ref_SBP, 'SBP Bland–Altman');

figure;
bland_altman(all_est_DBP, all_ref_DBP, 'DBP Bland–Altman');

save("BP_results.mat","Results");

fprintf("Pipeline complete! Saved to BP_results.mat\n");

%% =======================================================================
%   PART 2 — SAMPLE-LEVEL SUMMARY (30 samples)
% =======================================================================
total_rejected = 0;
for k = 1:numSamples
    SBPvals = [];
    DBPvals = [];
    REF_SBPvals = [];
    REF_DBPvals = [];

    for w = 1:numWindows
        SBPvals(end+1)    = Results(k).window(w).SBP_est;
        DBPvals(end+1)    = Results(k).window(w).DBP_est;
        REF_SBPvals(end+1)= Results(k).window(w).ref_SBP;
        REF_DBPvals(end+1)= Results(k).window(w).ref_DBP;
        
        total_rejected = total_rejected + Results(k).window(w).num_rej;

    end

    Results(k).sample.SBP_est_mean = mean(SBPvals);
    Results(k).sample.DBP_est_mean = mean(DBPvals);
    Results(k).sample.SBP_ref_mean = mean(REF_SBPvals);
    Results(k).sample.DBP_ref_mean = mean(REF_DBPvals);

    Results(k).sample.SBP_error = mean(SBPvals - REF_SBPvals);
    Results(k).sample.DBP_error = mean(DBPvals - REF_DBPvals);

    Results(k).sample.SBP_window_variance = var(SBPvals);
    Results(k).sample.DBP_window_variance = var(DBPvals);

end


%% =======================================================================
%   PART 3 — GLOBAL ACCURACY METRICS
% =======================================================================
valid = ~(isnan(all_est_SBP) | isnan(all_ref_SBP));

SBP_errors = all_est_SBP(valid) - all_ref_SBP(valid);
DBP_errors = all_est_DBP(valid) - all_ref_DBP(valid);

SBP_MAE = mean(abs(SBP_errors));
DBP_MAE = mean(abs(DBP_errors));

SBP_bias = mean(SBP_errors);
DBP_bias = mean(DBP_errors);

SBP_SD = std(SBP_errors);
DBP_SD = std(DBP_errors);



fprintf("\n=== GLOBAL ACCURACY SUMMARY ===\n");
fprintf("SBP  MAE = %.2f mmHg | Bias = %.2f | SD = %.2f\n", SBP_MAE, SBP_bias, SBP_SD);
fprintf("DBP  MAE = %.2f mmHg | Bias = %.2f | SD = %.2f\n", DBP_MAE, DBP_bias, DBP_SD);
fprintf("total rejected beats/pulse seg = %d  \n", total_rejected)
fprintf("total num detected beats = %d \n", num_pulses)
fprintf("total num detected segments = %d \n", num_segments)

%% =======================================================================
%   PART 4 — EXPORT RESULTS TO EXCEL
% =======================================================================
windowTable = [];
sampleTable = [];

for k = 1:numSamples
    for w = 1:numWindows
        row = Results(k).window(w);
        windowTable = [windowTable; ...
            table(k,w,row.SBP_est,row.ref_SBP,row.SBP_err,row.DBP_est,row.ref_DBP,row.DBP_err)];
    end

    sample = Results(k).sample;
    sampleTable = [sampleTable; ...
        table(k, sample.SBP_est_mean, sample.SBP_ref_mean, sample.SBP_error, ...
                 sample.DBP_est_mean, sample.DBP_ref_mean, sample.DBP_error)];
end

writetable(windowTable, "BP_window_results.xlsx");
writetable(sampleTable, "BP_sample_summary.xlsx");


%% =======================================================================
%   PART 5 — Improved Beat-to-Beat Visualization (SBP & DBP)
% =======================================================================

for k = 1:numSamples
    figure("Name",sprintf("Sample %d — Beat-to-Beat BP",k), "Color","w");

    all_SBP = []; all_DBP = []; all_t = [];
    ref_SBP = []; ref_DBP = []; ref_t = [];

    for w = 1:numWindows
        SBPvec = Results(k).window(w).SBP_beat;
        DBPvec = Results(k).window(w).DBP_beat;
        % Convert local 0–15s timestamps to global 0–60s timeline
        times_local = Results(k).window(w).peakTimes;
        t_offset = (w-1)*20 + 0;       % window start + 5 seconds skip % replace 0 with 5 if need be 
        times = times_local + t_offset;

        % Skip empty windows
        if isempty(SBPvec) || isempty(times), continue; end

        % Force equal lengths
        L = min([length(SBPvec), length(DBPvec), length(times)]);
        SBPvec = SBPvec(1:L);
        DBPvec = DBPvec(1:L);
        times  = times(1:L);

        all_SBP = [all_SBP; SBPvec(:)];
        all_DBP = [all_DBP; DBPvec(:)];
        all_t   = [all_t;   times(:)];

        % Reference BP for this window
        ref_SBP = [ref_SBP; Results(k).window(w).ref_SBP];
        ref_DBP = [ref_DBP; Results(k).window(w).ref_DBP];

        % Reference timestamps (flat within window)
        ref_t = [ref_t; repmat(times(1), length(SBPvec), 1)];
    end

   %% --- Plot SBP ---
    subplot(2,1,1);
    hold on;
    
    % Plot beat-to-beat SBP
    plot(all_t, all_SBP, 'b.-', "DisplayName","SBP estimated");
    
    % Plot reference SBP as flat lines per window
    for w = 1:numWindows
        t0 = (w-1)*20;     % window start
        t1 = w*20;         % window end
        rS = Results(k).window(w).ref_SBP;
        plot([t0 t1], [rS rS], 'r-', "LineWidth", 2, "DisplayName","Reference SBP");
    end
    
    title(sprintf("Sample %d — SBP Beat-to-Beat",k));
    xlabel("Time (s)");
    ylabel("SBP (mmHg)");
    grid on;
    
    % Window boundaries
    xline(20,'k--');
    xline(40,'k--');
    
    legend("Location","best");


    %% --- Plot DBP ---
    subplot(2,1,2);
    hold on;
    
    plot(all_t, all_DBP, 'g.-', "DisplayName","DBP estimated");
    
    for w = 1:numWindows
        t0 = (w-1)*20;
        t1 = w*20;
        rD = Results(k).window(w).ref_DBP;
        plot([t0 t1], [rD rD], 'r-', "LineWidth", 2, "DisplayName","Reference DBP");
    end
    
    title(sprintf("Sample %d — DBP Beat-to-Beat",k));
    xlabel("Time (s)");
    ylabel("DBP (mmHg)");
    grid on;
    
    xline(20,'k--');
    xline(40,'k--');
    
    legend("Location","best");

end

%% =======================================================================
%   PART 6 — Heatmap of Beat-Level Errors (Corrected)
% =======================================================================

for k = 1:numSamples
    
    figure("Name",sprintf("Sample %d — Beat-to-Beat SBP Error",k), ...
        "Color","w");
    hold on;

    all_err = [];
    all_t   = [];

    for w = 1:numWindows
        
        SBPvec = Results(k).window(w).SBP_beat;
        times_local = Results(k).window(w).peakTimes;
        ref_val = Results(k).window(w).ref_SBP;

        % Skip windows with no data
        if isempty(SBPvec) || isempty(times_local)
            continue;
        end

        % Ensure matching lengths
        L = min(length(SBPvec), length(times_local));
        SBPvec = SBPvec(1:L);
        times_local = times_local(1:L);

        % --- CRITICAL FIX: global timeline offset ---
        t_offset = (w-1)*20 + 0;   % skip first 5 seconds %replace 0 with 5 if needed
        times_global = times_local + t_offset;

        % --- Compute error for each beat ---
        err = SBPvec - ref_val;

        all_t = [all_t; times_global(:)];
        all_err = [all_err; err(:)];
    end

    % --- Safe plot check ---
    if isempty(all_t)
        warning("No valid beats for sample %d, skipping heatmap.", k);
        continue;
    end

    % --- Plot heatmap as colored scatter ---
    scatter(all_t, all_err, 45, all_err, 'filled');
    colorbar; colormap jet;

    title(sprintf("Sample %d — Beat-to-Beat SBP Error",k));
    xlabel("Time (s)");
    ylabel("SBP Error (mmHg)");

    % Window markers
    xline(20,'k--');
    xline(40,'k--');

    grid on;
end


%% =======================================================================
%   PART 7 — Beat Count Quality Plot
% =======================================================================

beatCount = zeros(numSamples, numWindows);

for k = 1:numSamples
    for w = 1:numWindows
        beatCount(k,w) = Results(k).window(w).stats.beat_count;
    end
end

figure("Name","Beat Count per Sample Window","Color","w");
imagesc(beatCount);
colorbar; colormap turbo;
xlabel("Window (1=0-20s, 2=20-40s, 3=40-60s)");
ylabel("Sample");
title("Beat Count Heatmap (Higher = Better Detection)");


%% =======================================================================
%   PART 8 — Boxplots of SBP/DBP Window Errors
% =======================================================================

SBP_errors = [];
DBP_errors = [];

for k = 1:numSamples
    for w = 1:numWindows
        SBP_errors(end+1) = Results(k).window(w).SBP_err;
        DBP_errors(end+1) = Results(k).window(w).DBP_err;
    end
end

figure("Color","w");
subplot(1,2,1); boxplot(SBP_errors);
title("SBP Window-Level Errors"); ylabel("Error (mmHg)");

subplot(1,2,2); boxplot(DBP_errors);
title("DBP Window-Level Errors"); ylabel("Error (mmHg)");


%% SAVE EVERYTHING
save("BP_results.mat","Results");
fprintf("\nAll results saved.\n");
