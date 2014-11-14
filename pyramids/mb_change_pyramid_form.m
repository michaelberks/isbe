%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_pyr pind] = mb_change_pyramid_form(pyramid, pind, pyr_type)

if nargin < 3
    pyr_type = 's';
end

switch pyr_type
    
    case 's'
        
%%%%%%%%%%%%%%%%%%%%
% Steerable pyramids
%%%%%%%%%%%%%%%%%%%%
if iscell(pyramid)
    %Pyrmaid is in Rose form - change to Simmoncelli
    %new pyramind is single column vector; could pre-allocate but can't be
    %arsed - we only do this once so performance isn't key;
    new_pyr = [];
    
    % Calculate num_levels and num_oris
    [num_levels num_oris] = size(pyramid);
    num_levels = num_levels - 2;
    
    %Check if we need to form pind
    pind = zeros(num_levels*num_oris +2, 2);
    count = 1;
    
    for lev = 1:num_levels+2
        
        for ori = 1:num_oris
            
            %levels comprised of num_oris sub-band orientations
            new_pyr = [new_pyr; pyramid{lev, ori}(:)]; %#ok
            
            %fill in pind if necessary
            pind(count, :) = size(pyramid{lev, ori});
            count = count+1;
            
            if lev == 1 || lev == num_levels+2
                % high/low band pass - no sub-band orientations
                break;
            end
        end
            
    end
        
    
else
    
    % Calculate num_levels and num_oris
    num_levels = length(unique(pind, 'rows')) - 1;
    num_oris = (size(pind, 1) - 2) / num_levels;
    
    %Pyramid is in Simmoncelli form - change to Rose
    new_pyr = cell(num_levels+2, num_oris);
    band_idx = cumsum(prod(pind, 2));

    for lev = 1:num_levels+2
        if lev == 1
            % high band pass - no sub-band orientations
            new_pyr{lev, 1} = reshape(pyramid(1:band_idx(1)), pind(1,1), pind(1,2));
        elseif lev == num_levels+2
            %low band pass - no sub-band orientations
            new_pyr{lev, 1} = reshape(pyramid(band_idx(end-1)+1:end), pind(end,1), pind(end,2));
        else
            %levels comprised of num_oris sub-band orientations
            for ori = 1:num_oris
                band = num_oris*(lev-2) + ori + 1;
                new_pyr{lev, ori} = reshape(pyramid(band_idx(band-1)+1:band_idx(band)),...
                    pind(band,1), pind(band,2));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%

    case 'g'

%%%%%%%%%%%%%%%%%%%%
% Gaussian pyramids
%%%%%%%%%%%%%%%%%%%%
if iscell(pyramid)
    %Pyrmaid is in Rose form - change to Simmoncelli
    %new pyramind is single column vector; could pre-allocate but can't be
    %arsed - we only do this once so performance isn't key;
    new_pyr = [];
    
    % Calculate num_levels and num_oris
    
    num_levels = size(pyramid, 1);
    
    %Check if we need to form pind
    pind = zeros(num_levels, 2);
    
    for lev = 1:num_levels
         
        %levels comprised of num_oris sub-band orientations
        new_pyr = [new_pyr; pyramid{lev}(:)]; %#ok

        %fill in pind if necessary
        pind(lev, :) = size(pyramid{lev});
            
    end
           
else
    
    % Calculate num_levels and num_oris
    num_levels = size(pind, 1);
    
    %Pyramid is in Simmoncelli form - change to Rose
    new_pyr = cell(num_levels, 1);

    for lev = 1:num_levels
        
        new_pyr{lev} = ...
            reshape(pyramid(1:prod(pind(lev,:))), pind(lev,1), pind(lev,2));
        pyramid(1:prod(pind(lev,:))) = [];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
end
