function [fg_indices, bg_indices] = get_precomputed_indices(images)

%If we've already computed which samples to take from this image load these
%now
if isfield(images, 'precomputed_indices')
    idx = load(images.precomputed_indices);
    fg_indices = idx.fg_indices;
    bg_indices = idx.bg_indices;
else
    fg_indices = [];
    bg_indices = [];
end


