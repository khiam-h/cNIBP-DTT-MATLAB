function stats = summarize_features(features)

    names = ["IBI","HR","Amp","RiseT","DecayT","UpSlope","DownSlope","Area", "DTT"];
    stats = struct();

    for i = 1:length(names)
        x = features(:,i);

        stats.(names(i)+"_mean")   = mean(x);
        stats.(names(i)+"_std")    = std(x);
        stats.(names(i)+"_median") = median(x);
        stats.(names(i)+"_IQR")    = iqr(x);
        stats.(names(i)+"_CV")     = std(x)/mean(x);
        stats.(names(i)+"_skew")   = skewness(x);
        stats.(names(i)+"_kurt")   = kurtosis(x);
    end

    stats.beat_count = size(features,1);
end

%ignore IBI/riseT/DecayT