function sw_lookup = sw_lookup_table(wedgevals)

% need a lookup table from pixel value (i_pix) to x_sw
sw_lookup = zeros(256,1);

% find number of steps
dim = size(wedgevals);
numsteps=dim(1);

% loop over i_pix i=1-256 corresponds to i_pix = 0-255
for i=1:256
    % if i_pix is off end of wedge then set x_sw to zero or 25/35mm
    if (i-1)>wedgevals(1)
        sw_lookup(i)=0;
    elseif (i-1)==wedgevals(1)
        sw_lookup(i)=1;
    elseif (i-1)<=wedgevals(numsteps)
        sw_lookup(i)=numsteps;    
    else
        i_step = 1;
        wedgevals(i_step)
        while wedgevals(i_step)>(i-1)
            i_step=i_step+1;
        end
        
        i
        i_step
        size(sw_lookup)
        size(wedgevals)
       
        sw_lookup(i) = ((i-1) - wedgevals(i_step-1))/(wedgevals(i_step)-wedgevals(i_step-1)) + (i_step-1);
    end
end

% move all x_sw vals which are above the highest we can measure back down to that highest value:
for i=1:256
    if(sw_lookup(i)<numsteps)
        max_val = sw_lookup(i);
        max_step = i;
        break;
    end
end

sw_lookup(1:i-1) = max_val;