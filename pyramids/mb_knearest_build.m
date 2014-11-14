function [knearest_data] = mb_knearest_build(varargin)
%
% MB_KNEAREST_BUILD Synthesise 
%
% MB_KNEAREST_BUILD uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Return Value:
%
%   MB_KNEAREST_BUILD returns KNEAREST_DATA.
%
% Notes: Although the base of this code provides general k-nearest
% functionality, it has been tweaked for my specific use with image
% pyramids.
%
% References:
%
% See also: MB_KNEAREST_SEARCH

% Unpack the arguments:
args = u_packargs(varargin, '0', ... % strict mode
		  {...  % The mandatory arguments
          'Image1',...
          'Image2',...
          },... % The optional arguments
          'BothLevels', false, ...
          'WindowSize1', 5,...
          'WindowSize2', 9,...
          'Method', 'standard',...
          'ROISubscripts', []...
          );

clear varargin;



%See if we've been given an ROI from which to select data
if isempty(args.ROISubscripts)
    % If not use as much of image as possible
    
    %get dimensions of smaller (coarser) image
    [rows1, cols1] = size(args.Image1);
    
    %get overlap required at edge of image
    overlap = (args.WindowSize1 + 1) / 2;
    
    [rows cols] = mesgrid(overlap:rows1 - overlap, overlap:cols1 - overlap);
    args.ROISubscripts = [rows(:) cols(:)];
     
end



%work out data dimensions so we can pre-allocate
if ~args.BothLevels
    data_width = args.WindowSize1^2;
    data_length = size(args.ROISubscripts, 1);
    knearest_data.Values = zeros(data_length, 4);
else
    data_width = args.WindowSize1^2 + args.WindowSize2^2;
    data_length = 4*size(args.ROISubscripts, 1);
    knearest_data.Values = zeros(data_length, 1);
end

knearest_data.Data = zeros(data_length, data_width);
knearest_data.Subscripts = zeros(data_length, 2);

idx = 1;
for ii = 1:size(args.ROISubscripts, 1)
    
    row = args.ROISubscripts(ii,1);
    col = args.ROISubscripts(ii,2);
    if ~args.BothLevels

        knearest_data.Data(idx, :) = reshape(...
            sample_window(args.Image1, args.WindowSize1, row, col), 1, []);

        knearest_data.Subscripts(idx, :) = [row col];
        knearest_data.Values(idx, 1) = args.Image2(2*row-1, 2*col-1);
        knearest_data.Values(idx, 2) = args.Image2(2*row-1, 2*col);
        knearest_data.Values(idx, 3) = args.Image2(2*row, 2*col-1);
        knearest_data.Values(idx, 4) = args.Image2(2*row, 2*col);
        idx = idx+1;
    else
        sample1 = reshape(sample_window(...
            args.Image1, args.WindowSize1, row, col), 1, []);

        %2i-1, 2j-1
        knearest_data.Data(idx, 1:args.WindowSize1^2) = sample1;

        knearest_data.Data(idx, args.WindowSize1^2+1:end) = ...
            reshape(sample_window(args.Image2, args.WindowSize2,...
                2*row-1, 2*col-1), 1, []);
        knearest_data.Subscripts(idx, :) = [2*row-1, 2*col-1];
        knearest_data.Values(idx, 1) = args.Image2(2*row-1, 2*col-1);
        idx = idx+1;

        %2i-1, 2j
        knearest_data.Data(idx, 1:args.WindowSize1^2) = sample1;

        knearest_data.Data(idx, args.WindowSize1^2+1:end) = ...
            reshape(sample_window(args.Image2, args.WindowSize2,...
                2*row-1, 2*col), 1, []);
        knearest_data.Subscripts(idx, :) = [2*row-1, 2*col];
        knearest_data.Values(idx, 1) = args.Image2(2*row-1, 2*col);
        idx = idx+1;

        %2i, 2j-1
        knearest_data.Data(idx, 1:args.WindowSize1^2) = sample1;

        knearest_data.Data(idx, args.WindowSize1^2+1:end) = ...
            reshape(sample_window(args.Image2, args.WindowSize2,...
                2*row, 2*col-1), 1, []);
        knearest_data.Subscripts(idx, :) = [2*row, 2*col-1];
        knearest_data.Values(idx, 1) = args.Image2(2*row, 2*col-1);
        idx = idx+1;

        %2i, 2j
        knearest_data.Data(idx, 1:args.WindowSize1^2) = sample1;

        knearest_data.Data(idx, args.WindowSize1^2+1:end) = ...
            reshape(sample_window(args.Image2, args.WindowSize2,...
                2*row, 2*col), 1, []);
        knearest_data.Subscripts(idx, :) = [2*row, 2*col];
        knearest_data.Values(idx, 1) = args.Image2(2*row, 2*col);
        idx = idx+1;

    end
    if (any(isnan(knearest_data.Data(:))));
        warning('NaN value found, check ROI boundaries'); %#ok
    end
end

if ~args.BothLevels
    knearest_data.Mean1 = mean(knearest_data.Data(:));
    knearest_data.Std1 = std(knearest_data.Data(:));
    knearest_data.Data = (knearest_data.Data - knearest_data.Mean1) / ...
        knearest_data.Std1;
else
    d1 = knearest_data.Data(:,1:args.WindowSize1^2);
    d2 = knearest_data.Data(:,args.WindowSize1^2 + 1:end);
    
    knearest_data.Mean1 = mean(d1(:));
    knearest_data.Std1 = std(d1(:));
    d1 = (d1 - knearest_data.Mean1) / ...
        knearest_data.Std1;
    
    knearest_data.Mean2 = mean(d2(:));
    knearest_data.Std2 = std(d2(:));
    d2 = (d2 - knearest_data.Mean2) / ...
        knearest_data.Std2;
    
    knearest_data.Data = [d1 d2];
end

switch args.Method
    case 'standard'
        %data is as we need it
        return
        
    case 'tree'
        %build tsvq tree from data
        [knearest_data.Data knearest_data.TreeIndices]=...
            mb_tsvq_build('Data', knearest_data.Data);
        
end
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        