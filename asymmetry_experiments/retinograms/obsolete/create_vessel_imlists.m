function imlists = create_vessel_imlists(args)
% Create a structure that holds the filenames for all images of interest to
% a vessel predictor

d = dir([args.image_dir '/*.mat']);
if ~isempty(d), [imlists(1:numel(d)).raw] = deal(d.name); end

d = dir([args.foveal_mask_dir '/*.mat']);
if ~isempty(d), [imlists.f_mask] = deal(d.name); end

d = dir([args.vessel_mask_dir '/*.mat']);
if ~isempty(d), [imlists.v_mask] = deal(d.name); end

d = dir([args.ori_dir '/*.mat']);
if ~isempty(d), [imlists.ori_gt] = deal(d.name); end

d = dir([args.width_dir '/*.mat']);
if ~isempty(d), [imlists.width] = deal(d.name); end
   

% Check which images are selected - we'll assume if images are selected the
% user has managed to index images within the corrcet range
if ~isempty(args.selected_images)
    imlists = imlists(args.selected_images);
end
    