function [frames_h5] = load_hdf5_frames(hdf5_file, frames_path, rgb, frames_idx, rotate, save_dir)
%LOAD_HDF5_FRAMES *Insert a one line summary here*
%   [frames_rgb] = load_hdf5_frames(hdf5_file, frame_dims, frame_idx)
%
% Inputs:
%      hdf5_file - *Insert description of input variable here*
%
%      frames_path - *Insert description of input variable here*
%
%      frames_idx - *Insert description of input variable here*
%
%
% Outputs:
%      frames_rgb - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ~exist('frames_path', 'var') || isempty(frames_path)
    frames_path = '/frames/data';
end

if ~exist('rgb', 'var') || isempty(rgb)
    rgb = false;
end

frames_h5 = h5read(hdf5_file, frames_path);

if ~exist('frames_idx', 'var') || isempty(frames_idx)
    frames_idx = 1:size(frames_h5,3);
end

if rgb
    frames_r = permute(frames_h5(3:3:end,:,frames_idx), [2 1 4 3]);
    frames_g = permute(frames_h5(2:3:end,:,frames_idx), [2 1 4 3]);
    frames_b = permute(frames_h5(1:3:end,:,frames_idx), [2 1 4 3]);
    frames_h5 = cat(3, frames_r, frames_g, frames_b);
else
    frames_h5 = permute(frames_h5, [2 1 3]);
end

if ~exist('rotate', 'var') || isempty(rotate)
    rotate = 2;
end
if rotate
    frames_h5 = rot90(frames_h5, rotate);
end

if ~exist('save_dir', 'var')
    save_dir = [];
end

if ~isempty(save_dir)
    create_folder(save_dir);
    for i_fr = 1:length(frames_idx)
        fr = frames_idx(i_fr);
        if rgb
            imwrite(uint8(frames_h5(:,:,:,fr)), [save_dir '\frame' zerostr(fr, 5) '.png']);
        else
            imwrite(uint8(frames_h5(:,:,fr)), [save_dir '\frame' zerostr(fr, 5) '.png']);
        end
    end
end
