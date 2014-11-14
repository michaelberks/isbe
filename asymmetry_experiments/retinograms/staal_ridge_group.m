function group_map = staal_ridge_group(vessel_prob, vessel_ori)

E_c = 3;
E_o = 0.95;
E_p = 0.9;
set_max = 25;
[rows cols] = size(vessel_prob);

if ~isreal(vessel_ori)
    %halve angles and normalise the vessel orientations
    vessel_theta = angle(vessel_ori)/2;
    
    %Get ridge pixels
    ridge_map = mb_non_maximal_supp(vessel_prob, vessel_theta) > 0;
else
    vessel_theta = vessel_ori;
    ridge_map = vessel_prob;
end

%Convert orientation to complex form
vessel_ori = exp(1i*vessel_theta);    

%initialise group map
group_map = zeros(size(vessel_prob));

%create patch to select pixels within desired radius of point
E_cc = ceil(E_c);
E_cv = -E_cc:E_cc;
xo = repmat(E_cv, 2*E_cc+1, 1);
yo = xo';
r = sqrt(xo.^2 + yo.^2);
xr = xo ./ r;
yr = -yo ./ r;
select_patch = r <= E_c;
select_patch(E_cc+1,E_cc+1) = 0;

go_on = true;
set_num = 0;
while go_on
    %Get list of potential ungrouped ridge pixels
    [yu xu] = find(~group_map & ridge_map);
    
    %If empty we're finished
    if isempty(yu);
        break
    end
    
    %Otherwise create a new set and pick a seed at random 
    set_num = set_num + 1;
    r_idx = ceil(length(yu)*rand);
    x_seed = xu(r_idx);
    y_seed = yu(r_idx);
    
    %Get  pixels grouped to seed
    [gx gy group_map] = get_groups(x_seed, y_seed, group_map, set_num);
    
    set_size = 1 + length(gx);
    
    %Greedily add group pixels
    while ~isempty(gx)
        
        %Take x,y points of first new grouped pixel
        xi = gx(1);
        yi = gy(1);
        
        %Remove grouped pixel from unchecked set
        gx(1) = [];
        gy(1) = [];
        
        %Check for new grouped pixels
        [gxi gyi group_map] = get_groups(xi, yi, group_map, set_num);
        
        %Add any to unchecked group
        gx = [gx; gxi]; %#ok
        gy = [gy; gyi]; %#ok
        
        %check we've not exceed the maximum number of pixels allowed per
        %set
        set_size = set_size + length(gxi);
        if set_size >= set_max
            break;
        end
    end
end

function [gxi gyi group_map] = get_groups(x_seed, y_seed, group_map, set_num)
    %First set seed points in group map and get orientation at seed
    group_map(y_seed, x_seed) = set_num;
    v_g = vessel_ori(y_seed, x_seed);
    
    %Get patch of local ugrouped ridge pixels within E_c
    ungroup_patch = ~group_map(y_seed+E_cv, x_seed+E_cv) & ...
        ridge_map(y_seed+E_cv, x_seed+E_cv) & select_patch;
    
    %if this is balnk return - there are no grouped pixels
    if ~any(ungroup_patch)
        gxi = [];
        gyi = [];
        return;
    end
    
    %Get orientation and r vectors of ungrouped pixels local region
    ori_patch = vessel_ori(y_seed+E_cv, x_seed+E_cv);
    v_u = ori_patch(ungroup_patch);
    xru = xr(ungroup_patch);
    yru = yr(ungroup_patch);
    gxi = x_seed + xo(ungroup_patch);
    gyi = y_seed + yo(ungroup_patch);
    
    %Check orientations match
    o_check = abs(real(v_u)*real(v_g) + imag(v_u)*imag(v_g)) >= E_o;
    
    %Check parallel
    p_check = abs(xru*real(v_g) + yru*imag(v_g)) >= E_p;
    
    %Grouped pixels pass both tests
    g_idx = o_check & p_check;
    
    %Return pixels that pass both tests
    gxi = gxi(g_idx);
    gyi = gyi(g_idx);
    
    %Add all these pixels to the group map
    group_map(sub2ind([rows cols], gyi, gxi)) = set_num;
    
end


end