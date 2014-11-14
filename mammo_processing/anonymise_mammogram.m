function anonymise_mammogram(mammo_dir, varargin)
%MB_MAKE_NEW_FUNCTION creates a new function with documentation template
%   [] = u_newfunction(function_name,ninput,i1,i2,...,in,o1,o2,...on)
%
%   inputs:
%      mammo_dir - String specifying directory containing mammograms to
%       anonymise
%   optional arguments pairs:
%
%      anonymised_dir - String specifying directory to save anonymised
%       mammograms to. If empty overwrites existing files
%
%      mammo_list - Structure of file names specifying mammograms to
%        anonymise. If empty will select all files of given 'image_type'
%
%      image_type - String specifying file extension to search for in input
%       directory if mammo_list not already specified
%
%   outputs:
%
%   
%   notes:

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 'anonymised_dir', [],...
             'mammo_list', [],...
             'image_type', '.bmp',...
             'anon_size', [350 85],...
             'right_def', 'R');

%Check the mammogram directory is filesep terminated         
if ~isempty(mammo_dir) && ~strcmp(mammo_dir(end), filesep)
    mammo_dir = [mammo_dir filesep];
end

%Get list of mammograms in directory if not already specified
if isempty(args.mammo_list)
    args.mammo_list = dir([mammo_dir, '*', args.image_type]);
end

%Set output directory
if isempty(args.anonymised_dir)
    args.anonymised_dir = mammo_dir;
end

%Load in each mammogram and anonymise patient label
for ii = 1:length(args.mammo_list);
    mammo = imread([mammo_dir args.mammo_list(ii).name]);
    [r c] = size(mammo);
    
    %check if left or right mammogram
    if ~isempty(strfind(args.mammo_list(ii).name, args.right_def))
        %right mammogram
        r1 = 1; r2 = args.anon_size(1);
        c1 = c - args.anon_size(2); c2 = c;
    else 
        %left mammogram
        r1 = r - args.anon_size(1); r2 = r;
        c1 = c - args.anon_size(2); c2 = c;
    end
    
    grey_val = mean(mean(mammo(r1:r2, c1:c2)));
    mammo(r1:r2, c1:c2) = grey_val;
    imwrite(mammo, [args.anonymised_dir args.mammo_list(ii).name]); 
    clear mammo;
    
end