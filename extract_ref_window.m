function ref_win = extract_ref_window(ref, t_start, t_end)

    idx = (ref.time >= t_start) & (ref.time < t_end);

    if any(idx)
        ref_win.SBP_mean = mean(ref.SBP(idx));
        ref_win.DBP_mean = mean(ref.DBP(idx));
        ref_win.HR_mean  = mean(ref.HR(idx));
        ref_win.valid = true;
    else
        ref_win.SBP_mean = NaN;
        ref_win.DBP_mean = NaN;
        ref_win.HR_mean = NaN;
        ref_win.valid = false;
    end
end
