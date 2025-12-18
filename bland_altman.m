function bland_altman(est, ref, titleStr)

    est = est(:);
    ref = ref(:);

    meanVals = (est + ref) / 2;
    diffVals = (est - ref);

    bias = mean(diffVals);
    sd   = std(diffVals);

    loa_upper = bias + 1.96*sd;
    loa_lower = bias - 1.96*sd;

    scatter(meanVals, diffVals, 40, 'filled'); hold on;
    yline(bias, 'r--', sprintf('Bias = %.2f', bias));
    yline(loa_upper, 'k--', sprintf('Upper LoA = %.2f', loa_upper));
    yline(loa_lower, 'k--', sprintf('Lower LoA = %.2f', loa_lower));

    xlabel('Mean of Estimated and Reference (mmHg)');
    ylabel('Difference (Estimate - Reference) (mmHg)');
    title(titleStr);
    grid on;
end
