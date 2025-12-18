function idx_valid = beat_filter(features, stats)

IBI      = features(:,1);
HR       = features(:,2);
Amp      = features(:,3);
RiseT    = features(:,4);
DecayT   = features(:,5);
UpSlope  = features(:,6);
DownSlope= features(:,7);
Area     = features(:,8);
DTT      = features(:,9);

idx_valid = true(size(IBI));

% Rule 1: Unphysiological IBI
idx_valid = idx_valid & (IBI > 0.4 & IBI < 1.7);

% Rule 2: Physiologically valid relaxation time
%idx_valid = idx_valid & (DecayT > 0);

% Rule 3: Rise time plausibility
%idx_valid = idx_valid & (RiseT > 0.05 & RiseT < 0.40);

% Rule 4: DownSlope plausibility (decay rate)
%idx_valid = idx_valid & (DownSlope > -2.0 & DownSlope < -0.01);

% Rule 5: Reasonable pulse amplitude
%idx_valid = idx_valid & (Amp > 0.02 & Amp < 5);

% Rule 6: Reasonable DTT
%DTT_uplimit = stats.DTT_median + 1*stats.DTT_std;
%DTT_lowlimit = stats.DTT_median - 1*stats.DTT_std;
%DTT  = idx_valid & (DTT > DTT_lowlimit & DTT < DTT_uplimit);

end
