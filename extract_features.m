function [features, peakTimes, stats] = extract_features(annot)

    sig = annot.resampSignal;
    t   = annot.resampTime;
    peaks = annot.sysPeakIndex;
    troughs = annot.diaTroughIndex;

    fs = 200;  % resampled by BP_annotate

    features = [];
    peakTimes = t(peaks);  % timestamps of each detected beat

    if length(peaks) < 3
        return;
    end

    for i = 1:length(peaks)-1
        fprintf("extract features run %d ... \n", (i) )
        if (peaks(i) == peaks(i+1)) || (troughs(i) == troughs(i+1))
            continue
        else 
            p1 = peaks(i);
            p2 = peaks(i+1);
            T2 = troughs(i+1);
    
            beat = sig(p1:p2);
    
            IBI = (t(p2) - t(p1));  % seconds
            HR  = 60 / IBI;
    
            amp = sig(p1) - min(beat);
    
            [~, maxLoc] = max(beat);
            riseT = maxLoc / fs;
            decayT = (length(beat) - maxLoc) / fs;
    
            d = diff(beat);
            upSlope   = max(d);
            downSlope = min(d);
    
            areaVal = trapz(beat);
            %fprintf("areaval = %d \n", areaVal)
    
            DTT = (t(T2) - t(p1));
            %fprintf("DTT = %d \n", DTT)
            
            test_array = [IBI HR amp riseT decayT upSlope downSlope areaVal DTT];
            features(end+1,:) = [IBI HR amp riseT decayT upSlope downSlope areaVal DTT];
        end
    end
    stats = summarize_features(features);
end
