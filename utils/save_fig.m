function [] = save_fig(fig_dir, fig_names, fig_ext, num_figs, close_figs)
%SAVE_FIG export a set of open figures, uses exportfig
%   [] = save_fig(fig_dir, fig_names, num_figs)
%
% Inputs:
%      fig_dir - directory into which figures are saved
%
%      fig_names - either a cell string with a name for each figure, or a
%      string used as the figure root, to which the figure number is
%      appended. Eg if fig_names = 'fig', the figures will be save as
%      'fig001.png', 'fig002.png' etc. Supplying a cell string shorter than
%      the number of figures to save will cause an error when it tries to
%      save nth+1 image
%
%      fig_ext - extension (including .) of used to specify output type.
%      Default '.png'
%
%      num_figs - Number of figures to save, note the ordering of figures
%      may not be guaranteed, I suspect they come in order of when they
%      opended but should probably check. If not supplied or a number more 
%      than theall open figs
%      will
%
%      close_figs - flag to choose whether to close each figure after
%      saving. Default - false
%
%
% Outputs: none
%
% Example: [] = save_fig('C:\temp_figs', 'my_fig')
%
% Notes:
%
% See also: EXPORTFIG
%
% Created: 12-Oct-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
fig_handles =  findobj('type','figure');

if ~exist('close_figs', 'var')
    close_figs = false;
end

if exist('num_figs', 'var')
    num_figs = min(length(fig_handles), num_figs);
else
    num_figs = length(fig_handles);
end

if ~exist('fig_ext', 'var') || isempty(fig_ext)
    fig_ext = '.png';
end

fig_dir = create_folder(fig_dir);

for i_fig = 1:num_figs
    if iscellstr(fig_names) 
        if length(fig_names) > i_fig
            fig_name = fullfile(fig_dir, [fig_names{i_fig} fig_ext]);
        else
            error(['Length of names list (' num2str(length(fig_names)) ...
                ') is smaller than number of figs to save (' num2str(num_figs) ')']);
        end
    elseif ischar(fig_names)
        fig_name = fullfile(fig_dir, [fig_names sprintf('%03d', i_fig) fig_ext]);
    else
        error(['fig_names must be a string or cell of strings']);
    end
    figure(fig_handles(i_fig));
    exportfig(fig_name);
    display(['Saved figure to ' fig_name]);
    
    if close_figs
        close(fig_handles(i_fig));
    end
end

