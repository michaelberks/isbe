function [] = mb_sample_valid_region(varargin)
%MB_SAMPLE_VALID_REGION *Insert a one line summary here*
%   [] = mb_sample_valid_region(varargin)
%
% MB_SAMPLE_VALID_REGION uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Apr-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0',...
        {...
        'WindowSize', ...
        'MammoDir', ...
        'SaveDir' ...
        },...
        'NameStem', [],...
        'NumRegions', 1,...
        'MammoFiles', [], ...
        'MaskDir', [],...
        'MaskFiles',[],...
        'ResizeFactor', 4, ...
        'ImageFormat', '.mat' ...
        );

%Check mammo directory is filesep terminated
if ~strcmp(args.MammoDir(end), filesep)
    args.MammoDir = [args.MammoDir filesep];
end

%Get mammo filenames if not specified
if isempty(args.MammoFiles)
    args.MammoFiles = dir([args.MammoDir, '*', args.ImageFormat]);
end

%Make sure we have a valid directory for masks
if isempty(args.MaskDir)
    args.MaskDir = args.MammoDir;
elseif ~strcmp(args.MaskDir(end), filesep)
    args.MaskDir = [args.MaskDir filesep];
end

%Get mask filenames if not specified
if isempty(args.MaskFiles)
    args.MaskFiles = dir([args.MaskDir, '*mask*']);
end

%Check file lists are same length
if length(args.MaskFiles) ~= length(args.MammoFiles)
    error('The number of masks and mammograms did not match');
end

%Make sure we have a valid directory to save samples to
if isempty(args.SaveDir)
    args.SaveDir = args.MammoDir;
elseif ~strcmp(args.SaveDir(end), filesep)
    args.SaveDir = [args.SaveDir filesep];
end

%If we going to save the sampled regions using a name stem and sequential
%numbering we need to set up a counter here
if ~isempty(args.NameStem)
    curr_im = 1;
end

for ii = 1:length(args.MammoFiles)
    
    %load the mask, reduce size and clear original
    mask = u_load([args.MaskDir, args.MaskFiles(ii).name]);
    mask_small = imresize(mask, 1/args.ResizeFactor);
    clear mask;
    
    %compute reduced window size and erode small mask
    mask_small([1 end], :) = 0;
    mask_small(:, [1 end]) = 0;
    window_size_small = ceil(args.WindowSize / args.ResizeFactor) + 1;
    mask_eroded = imerode(mask_small, ones(window_size_small));
    clear mask_small;
    
    %get valid rows and cols
    [rows cols] = find(mask_eroded);
    clear mask_eroded;
    
    %Check there are valid positions to sample from, if not, continue to
    %next image
    if isempty(rows)
        continue;
    end
    
    %load mammogram
    if strcmp(args.ImageFormat, '.mat')
        mammo = u_load([args.MammoDir, args.MammoFiles(ii).name]);
    else
        mammo = imread([args.MammoDir, args.MammoFiles(ii).name]);
    end
   
    
    for jj = 1:args.NumRegions
        
        %Randomly select a valid location
        idx = ceil(rand * length(rows));
        row = rows(idx); col = cols(idx);

        %expand row and col to full size
        row = round(row * args.ResizeFactor) + 1;
        col = round(col * args.ResizeFactor) + 1;

        %Extract region
        sampled_window = sample_window(mammo, args.WindowSize, row, col); %#ok
        
        %Set up the save name and save
        if ~isempty(args.NameStem)
            save_name = [args.NameStem zerostr(curr_im, 5) '.mat'];
            curr_im = curr_im + 1;
        else
            %Name region with the mammo name and information on how it was
            %sampled
            save_name = [args.MammoFiles(ii).name(1:end-4),...
                '_', num2str(jj),...
                '_', num2str(args.WindowSize), ...
                '_', num2str(row), ...
                '_', num2str(col) '.mat'];
        end    
        save([args.SaveDir, save_name], 'sampled_window');
    end
    clear mammo;
end
