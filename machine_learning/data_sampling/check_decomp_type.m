function decomp_type = check_decomp_type(decomp_type)

if ~nargin && ~nargout
    test_func;
    return;
end

%This function actually does some work now - first we make sure we have
%cell string
if ~iscell(decomp_type)
    decomp_type = cellstr(decomp_type);
end

%Now we check each decomp type in the cell is valid
n_types = length(decomp_type);
for ii = 1:n_types
    switch decomp_type{ii}
        case {'pixel',... % dumbest
              'linop', ... % template matching
              'g', 'g1d', 'g2d', 'g2di', 'g2da', 'g2dia', 'h2d', 'h2di', 'h2da', 'h2dia', 'haar', ... % derivs
              'dt', 'mono', 'gabor', 'gabori'... % with phase
             }
          % do nothing
        case 'clover'
            error('Clover is no longer a supported decomposition type, please use g2d instead');
            
        %Do conversions of old combination types - we can assume in this
        %case it's just a single string (otherwise it's not safe to use anyway)        
        case 'g12d'
            decomp_type{1} = 'g1d';
            decomp_type{2} = 'g2d';
            warning('ASYM:decomp_type', 'Converting old decomposition format');
            
        case 'dtg2'
            decomp_type{1} = 'dt';
            decomp_type{2} = 'g2d';
            warning('ASYM:decomp_type', 'Converting old decomposition format');
            
         case 'g2dh'
            decomp_type{1} = 'g2d';
            decomp_type{2} = 'h2d';
            warning('ASYM:decomp_type', 'Converting old decomposition format');
            
         case 'g2dg'
            decomp_type{1} = 'g';
            decomp_type{2} = 'g2d';
            warning('ASYM:decomp_type', 'Converting old decomposition format');
            
         case 'gaborg'
            decomp_type{1} = 'gabor';
            decomp_type{2} = 'g';
            warning('ASYM:decomp_type', 'Converting old decomposition format');

        otherwise
            error(['Unknown decomposition type: ',decomp_type{ii}]);
    end
end

function test_func

%All the old valid types before we swapped to the new cellstr format
for decomp_type = {'pixel',... % dumbest
          'linop', ... % template matching
          'g1d', 'g2d', 'g2di', 'g2da', 'haar', 'g2dg', ... % derivs
          'dt', 'mono', 'g12d', 'g2dh', 'gabor', 'gaborg',... % with phase
          'dtg2' % suggested by dumb reviewer
         }
     display(check_decomp_type(decomp_type));
end

