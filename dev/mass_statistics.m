function [S] = mass_statistics(varargin)

%
% MB_MASS_STATISTICS Compute summary statistics for a set of subtracted
% mammographic masses
%
%
% MB_MASS_STATISTICS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%   'MassNames' - structure of filenames to a set of masses
%   'MassPath' - Path to the MassNames directory 
%   'Masses' - structure of masses (exactly one of Masses or MassNames must
%     be set
%   'VarWindowSize' - size of window in which to calculate local
%     texture variance
%   'Statistics' - user specified list of statistics to return (not
%     implemented yet)
%   'SaveFile' - If non-empty, the pathname to save the returned structured
%
% Return Value:
%
%   MB_CLS_SELECTION returns S a structure containing the
%   following fields:
%    'Area' 
%      - Vector of area for each mass
%    'Perimeter' 
%      - Vector of perimeter for each masst
%    'TextureMean' 
%      - Vector of mean texture intensity for each mass
%     'TextureMin' 
%      - Vector of minimum texture intensity for each mass
%     'TextureMax' 
%      - Vector of maximum texture intensity for each mass
%     'TextureVar' 
%      - Vector of mean local variance of texture intensity for each mass

args = u_packargs(varargin, '0', ... % strict mode
		  {},... 
          'MassNames', [],... % The optional arguments
          'MassPath', [],...
          'Masses', [],...,
          'VarWindowSize', 11,...
          'SaveFile', []);

if ~isempty(args.Masses)
    N = length(args.Masses);
    load_masses = 0;
elseif ~isempty(args.MassNames)
    N = length(args.MassNames);
    load_masses = 1;
end

S.Area = zeros(N,1);
S.Perimeter = zeros(N,1);
S.TextureMean = zeros(N, 1);
S.TextureMin = zeros(N, 1);
S.TextureMax = zeros(N, 1);
S.TextureVar = zeros(N, 1);

for ii = 1:N
    if load_masses
        mass = u_load([args.MassPath, args.MassNames(ii).name]);
    else
        mass = args.Masses(ii);
    end
    
    mass_bw = roipoly(mass.subtract_ROI, mass.mass_outline(:,1),...
        mass.mass_outline(:,2));

    mass_rp = regionprops(bwlabel(mass_bw,4),...
        'Area', 'Centroid', 'Perimeter');
    S.Area(ii) = sum(mass_bw(:));
    S.Perimeter(ii) = mass_rp.Perimeter;

    [im_var im_mean var_map] = image_stats(mass.subtract_ROI, 11);

    S.TextureMean(ii) = mean(mass.subtract_ROI(mass_bw));
    S.TextureMin(ii) = min(mass.subtract_ROI(mass_bw));
    S.TextureMax(ii) = max(mass.subtract_ROI(mass_bw));

    S.TextureVar(ii) = im_var;
    S.TextureVar(ii) = mean(var_map(mass_bw));
%         
    clear mass;
end

%Main function end