function [hist_struc, return_val] = hfill1s (hist_struc, x, weight ) 	

%See if x is in valid range
if (x < hist_struc.xmax && x >= hist_struc.xmin )
    
    %Find nearest bin
    i = floor((x - hist_struc.xmin)/ hist_struc.xincr + 0.5);
    if (weight == 0.0) 
        return_val = hist_struc.counts(i);
        return;
    end

    %Update summary statistics
    hist_struc.sum = hist_struc.sum + weight;
    hist_struc.entries = hist_struc.entries + 1;

    hist_struc.mean = hist_struc.mean +...
        weight*(x-hist_struc.mean)/hist_struc.sum;
    hist_struc.mean2 = hist_struc.mean2 +...
        weight*(x*x-hist_struc.mean2)/hist_struc.sum;

    % correct numerical instability on first time through
    if (hist_struc.mean*hist_struc.mean > hist_struc.mean2)
        hist_struc.mean2 = hist_struc.mean*hist_struc.mean;
    end
    
    %Compute scaled distance of x from bin centre - will vary [-0.5 0.5]
    fac = (x-i*hist_struc.xincr-hist_struc.xmin) / hist_struc.xincr;
    if (i>0) 
        hist_struc.counts(i-1) = hist_struc.counts(i-1) + ...
            weight*(0.5-fac); 
    end
    if (i < hist_struc.xbins) 
        hist_struc.counts(i) = hist_struc.counts(i) +...
            weight*(0.5+fac);
    end
    return_val = (hist_struc.counts(i));
else

    if (x>hist_struc.xmax) 
        hist_struc.over = hist_struc.over + weight;
    else
        hist_struc.under = hist_struc.under + weight;        
    end
    return_val = 0.0;
end