function [outlierLCC, outlierRCC, outlierLML, outlierRML] = ...
    volpara_outlier_rejection(densityLCC, densityRCC, densityLML, densityRML)

minSsdToMaxSsdThreshold = 0.4; %Min SSD to max SSD ratio threshold
medianDensitySF = 3; %Median density scale factor
medianDensityTO = 0; %Median density threshold offset

if (densityLCC > 0.1 && densityRCC > 0.1 && densityLML > 0.1 && densityRML > 0.1)
    % IWW Set Density
    densityMedian = median([densityLCC, densityRCC, densityLML, densityRML]);
    ssdLCC = 0.67 * (sum([densityRCC, densityLML, densityRML].^2) - densityRCC * densityLML - densityRCC * densityRML - densityLML * densityRML);
    ssdRCC = 0.67 * (sum([densityLCC, densityLML, densityRML].^2) - densityLCC * densityLML - densityLCC * densityRML - densityLML * densityRML);
    ssdLML = 0.67 * (sum([densityRCC, densityLCC, densityRML].^2) - densityRCC * densityLCC - densityRCC * densityRML - densityLCC * densityRML);
    ssdRML = 0.67 * (sum([densityRCC, densityLML, densityLCC].^2) - densityRCC * densityLML - densityRCC * densityLCC - densityLML * densityLCC);

    minSsd = min([ssdLCC, ssdRCC, ssdLML, ssdRML]);
    maxSsd = max([ssdLCC, ssdRCC, ssdLML, ssdRML]);

    thresholdSsd = densityMedian * medianDensitySF + medianDensityTO;

    if maxSsd == 0
        minMaxVal = 100;
    else
       minMaxVal = (minSsd / maxSsd) ^ 0.5;
    end

    outlierFarEnough = false;
    if maxSsd > thresholdSsd 
        outlierFarEnough = true;
    end

    minMaxSpread = false;
    if minMaxVal < minSsdToMaxSsdThreshold 
        minMaxSpread = true;
    end

    %Reject if Yes and Yes
    outlierLCC = false; outlierRCC = false; outlierLML = false; outlierRML = false; 
    if (outlierFarEnough && minMaxSpread)
        if minSsd == ssdLCC, outlierLCC = true; end
        if minSsd == ssdRCC, outlierRCC = true; end
        if minSsd == ssdLML, outlierLML = true; end
        if minSsd == ssdRML, outlierRML = true; end
    end
end