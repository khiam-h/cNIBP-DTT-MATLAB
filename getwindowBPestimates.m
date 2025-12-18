function [SBP_est, DBP_est, BP_mmhg, timestamps, rejection,num_seg] = getwindowBPestimates(sig_win, annot,calibration_params, fs, stats)

    t = 1/fs;

    SBP_0 = calibration_params.cuff_ref.initial_SBP;
    DBP_0 = calibration_params.cuff_ref.initial_DBP;
    m_0 = calibration_params.m_0;
    c_0 = calibration_params.c_0;
    avg_eDBP = 0;
    
    systolic_peaks = annot.sysPeakIndex;
    diastolic_troughs = annot.diaTroughIndex;
    timestamps = annot.resampTime;
    filtered_signal_dB = annot.resampSignal;


    %segment pulse annotations within window
    pulse_segments = struct();

    %DTT limits
    DTT_med = stats.DTT_median;
    DTT_uplimit = stats.DTT_median + 2*stats.DTT_std;
    DTT_lowlimit = stats.DTT_median - 2*stats.DTT_std;
    IBI = 60/calibration_params.cuff_ref.initial_HR;

    success_pulse = 0;
    rejection = 0;
    num_seg = 0;
    DTT_array = [];
    
    for i = 1:(length(diastolic_troughs)-1)
        start_id = diastolic_troughs(i);
        end_id = diastolic_troughs(i+1);
        DTT = timestamps(end_id)-timestamps(systolic_peaks(i));
        cycle_len = timestamps(end_id)-timestamps(start_id);
        DTT_array(end+1) = DTT;
        
        if DTT > 0.4*IBI && DTT < 0.9*IBI
            %success_pulse = success_pulse+1;
            num_seg = success_pulse;
            success_pulse = success_pulse +1;
            if cycle_len > 0 %0.4 && cycle_len < 1.7
                pulse_segments(success_pulse).DTT = DTT;
                pulse_segments(success_pulse).start_id = start_id;
                pulse_segments(success_pulse).end_id = end_id;
                pulse_segments(success_pulse).DBP_s = filtered_signal_dB(end_id);
                pulse_segments(success_pulse).SBP_s = filtered_signal_dB(systolic_peaks(i));
                pulse_segments(success_pulse).waveform = filtered_signal_dB(start_id:end_id);
                pulse_segments(success_pulse).SBP_id = systolic_peaks(i);
                pulse_segments(success_pulse).DBP_id = diastolic_troughs(i);

            else
                pulse_segments(success_pulse).DTT = calibration_params.DTT_avg;
                pulse_segments(success_pulse).start_id = start_id;
                pulse_segments(success_pulse).end_id = end_id;
                pulse_segments(success_pulse).DBP_s = filtered_signal_dB(end_id);
                pulse_segments(success_pulse).SBP_s = filtered_signal_dB(systolic_peaks(i));
                pulse_segments(success_pulse).waveform = filtered_signal_dB(start_id:end_id);
                pulse_segments(success_pulse).SBP_id = systolic_peaks(i);
                pulse_segments(success_pulse).DBP_id = diastolic_troughs(i);
            end
            
        else
            fprintf("DTT for pulse %d is %d ... \n",i,DTT)
            rejection = rejection + 1;

        end
    end

    for i = (1:length(pulse_segments))
        DTT_limit = stats.DTT_mean + stats.DTT_std;
        

        DTT_t = pulse_segments(i).DTT;
        c_t = max(gradient(pulse_segments(i).waveform));
        if c_t > 0 & ~isempty(DTT_t)
            estimated_DBP = SBP_0 - (m_0*DTT_t*(c_0/c_t));   
        else
            %continue
            estimated_DBP = DBP_0;
        end

        pulse_segments(i).eDBP = estimated_DBP;
        

    end


    %% Map window
    scaled_win = annot.resampSignal *calibration_params.scale_factor;
    %get eDBP interpolated data
    intpl_start = pulse_segments(1).end_id; % start point of interpolation (first DBP)
    intpl_end = pulse_segments(end).end_id; % Last DBP
    span = timestamps(intpl_start):t:timestamps(intpl_end);

    eDBP_timestamp = timestamps(1,[pulse_segments.end_id]);
    eDBP_table = [pulse_segments.eDBP];

    rem_eDBP = 0;
    eDBP_table(eDBP_table == rem_eDBP) =[];
    xq = span;
    C = unique(eDBP_timestamp);

    % Check if the number of elements in C is less than in A
    if length(C) < length(eDBP_timestamp)
        disp('There are duplicate sample points.');
        disp(eDBP_timestamp)
    else
        disp('All sample points are unique.');
        disp(length(eDBP_table))
        disp(length(eDBP_timestamp))
    end

    interpolated_eDBP = interp1(eDBP_timestamp, eDBP_table,xq, "cubic");


    %obtain baseline shift of eDBP 
    baseline_estimate_eDBP = medfilt1(interpolated_eDBP, 2000 ,"omitnan","truncate");
    padsize_post = length(baseline_estimate_eDBP)+ (length(scaled_win)-pulse_segments(end).end_id);
    baseline_estimate_eDBP = paddata(baseline_estimate_eDBP, padsize_post,"Dimension","auto","FillValue",DBP_0,"Side","trailing");
    padsize_pre = length(baseline_estimate_eDBP)+ pulse_segments(1).end_id - 1;
    baseline_estimate_eDBP = paddata(baseline_estimate_eDBP, padsize_pre,"Dimension","auto","FillValue",DBP_0,"Side","leading");

    %obtain baseline of scaled signal DBP
    scaled_win_DBP_s = [pulse_segments.DBP_s] * calibration_params.scale_factor;
    interpolated_scaled_DBP_s = interp1(eDBP_timestamp,scaled_win_DBP_s,xq, "cubic");
    baseline_scaled_DBP_s = medfilt1(interpolated_scaled_DBP_s, 10,"omitnan","truncate");
    scaledsig_padsize_post = length(baseline_scaled_DBP_s)+ (length(scaled_win)-pulse_segments(end).end_id);
    baseline_scaled_DBP_s = paddata(baseline_scaled_DBP_s, scaledsig_padsize_post,"Dimension","auto","Pattern","edge","Side","trailing");
    scaledsig_padsize_pre = length(baseline_scaled_DBP_s)+ pulse_segments(1).end_id - 1;
    baseline_scaled_DBP_s = paddata(baseline_scaled_DBP_s, scaledsig_padsize_pre,"Dimension","auto","Pattern","edge","Side","leading");

    %baseline correction of signal
    baseline_diff = baseline_scaled_DBP_s - baseline_estimate_eDBP;

    if length(scaled_win) ~= length(baseline_diff)
        baseline_diff(end) = [];
        baseline_estimate_eDBP(end) = [];
        corrected_signal = scaled_win - baseline_diff;
    else
        corrected_signal = scaled_win - baseline_diff;
    end
    

    [eqv_trough, eSBP] = BP_annotate(corrected_signal,fs,0);


    SBP_est  = [];
    DBP_est = [];
    rej_pulse = [];
    
    
    for i = 1:length(eSBP)
        SBP_val = corrected_signal(eSBP(i));
        DBP_val = baseline_estimate_eDBP(eqv_trough(i));
        if (SBP_val < 90 || SBP_val > 180) || (DBP_val < 60 || DBP_val > 120) % orig values 60, 200, 40, 150
            continue
            
        else
            SBP_est(end+1) = SBP_val;
            DBP_est (end+1) = DBP_val;
        end
    end


    BP_mmhg = corrected_signal;
    bp_timestamps = timestamps;
end


