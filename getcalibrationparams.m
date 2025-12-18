function calibration_params = getcalibrationparams(signal, ref, fs)
    
    [diastolic_troughs,systolic_peaks,notchindex,dicroticindex, timestamps,filtered_signal_dB] = BP_annotate(signal,fs,0);

    % get pulse segments
    pulse_segments = struct();
    min_IBI = 60/300; % physiological limits max hr = 300
    ref_hr = ref.HR(1);
    IBI_ref = 60/ref_hr;
    n = 500;
    success = 0;
    for i = 1:n
        start_id = diastolic_troughs(i);
        end_id = diastolic_troughs(i+1);
        DTT = timestamps(end_id)-timestamps(systolic_peaks(i));

        SBP_i = filtered_signal_dB(systolic_peaks(i));
        DBP_i = filtered_signal_dB(diastolic_troughs(i+1));
        DBP_i1 = filtered_signal_dB(diastolic_troughs(i+2));

        cycle_len = timestamps(end_id)-timestamps(start_id);
        eHR = 60/cycle_len;
        fprintf("looking at pulse %d \n", i)
        if DTT > 0.4*IBI_ref && DTT < 0.8*IBI_ref 
            fprintf("adding pulse")
            disp(DTT)
            success = success + 1;
            pulse_segments(success).DTT = DTT;
            pulse_segments(success).start_id = start_id;
            pulse_segments(success).end_id = end_id;
            pulse_segments(success).DBP_s = filtered_signal_dB(end_id);
            pulse_segments(success).SBP_s = filtered_signal_dB(systolic_peaks(i));
            pulse_segments(success).waveform = filtered_signal_dB(start_id:end_id);
        else
            continue
        end
        
        if success == 6
            break
        end
    end    

    %obtain calibration params
    cuff_ref = struct();
    cuff_ref.initial_SBP = ref.SBP(1);
    cuff_ref.initial_DBP = ref.DBP(1);
    cuff_ref.initial_HR = ref.HR(1);
    calibration_params = struct();
    contractility = 0;
    m_0 = 0;
    PP_0 = cuff_ref.initial_SBP - cuff_ref.initial_DBP;
    scaling_factor = 0;
    PP_0_norm = PP_0/cuff_ref.initial_DBP;
    DTT_avg = 0;

    ppsi_check = [];
    pps_next = [];

    for i = (1:5)
        %get m_0
        fprintf("run %d \n", i )
        PP_s_i = (pulse_segments(i).SBP_s - pulse_segments(i).DBP_s);
        fprintf("run end \n")
        PP_s_next = pulse_segments(i).SBP_s -pulse_segments(i+1).DBP_s;
        DTT_i = pulse_segments(i).DTT;
        m_0 = m_0 + (PP_0/PP_s_i)*(PP_s_next/DTT_i);

        %get contractility
        derivative = gradient(pulse_segments(i).waveform);
        contractility = contractility + max(derivative);
        
        %scaling factor
        scaling_factor = scaling_factor + PP_s_i; 

        %DTT_avg
        DTT_avg = DTT_avg + DTT_i;

        ppsi_check(end+1) = PP_s_i;
        pps_next(end+1) = PP_s_next;
    end
    
    m_0 = abs(m_0/5);
    contractility = contractility/5;
    scaling_factor = PP_0/(scaling_factor/5);
    DTT_avg = DTT_avg/5;
    calibration_params.m_0 = m_0;
    calibration_params.c_0 = contractility;
    calibration_params.scale_factor = scaling_factor;
    calibration_params.DTT_avg = DTT_avg;
    calibration_params.cuff_ref = cuff_ref;
    calibration_params.pulse_seg = pulse_segments;
    calibration_params.ppsi = ppsi_check;
    calibration_params.ppnext = pps_next;

    fprintf("Calibration Parameters Saved ... \n")
end
