function [rr, cc, num_samples] = get_indices_from_winsize(rows, cols, win_size)

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[rr, cc, num_samples] = func(rows, cols, win_size);


%% Current version of function
function [rr, cc, num_samples] = func(rows, cols, win_size)

% rows and cols must be column vectors
rows = rows(:);
cols = cols(:);

if (win_size == 1)
    rr = rows;
    cc = cols;
end

%Make copies of sample rows and cols at positions of local window patch
half_win_size = floor(win_size / 2);
centre_offset = -half_win_size:half_win_size;

num_samples = length(rows);
centre_offset_mat = ones(num_samples,1)*centre_offset;
ones_vec = ones(1,win_size);

% Repeat column vectors and apply the offsets
rr = rows(:,ones_vec) + centre_offset_mat;
cc = cols(:,ones_vec) + centre_offset_mat;

% Replicate columns: [a,b,c] -> [a,b,c,a,b,c,a,b,c,...]
rr = kron(ones_vec, rr);

% Replicate columns: [a,b,c] -> [a,a,a,...,b,b,b,...,c,c,c]
cc = kron(cc, ones_vec);


%% Old code (for comparison during updates)
function [rr, cc, num_samples] = old_func(rows, cols, win_size)
%Make copies of sample rows and cols at positions of local window patch
pad_w = floor(win_size/2);
win_idx = -pad_w:pad_w;
num_samples = length(rows);
rr = repmat(rows*ones(1,win_size) + ones(num_samples,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_samples,1)*win_idx, ones(1,win_size));


%% Test script
function test_script()
clc;

% Create dummy inputs here
inds = ceil(rand(10,1)*10000);
[rows,cols] = ind2sub([100,100], inds);
win_size = 3;

% Get outputs from both old and new version of code
tic;
[rr1, cc1, num_samples1] = old_func(rows, cols, win_size);
t1 = toc;
tic;
[rr2, cc2, num_samples2] = func(rows, cols, win_size);
t2 = toc;

% Success true if outputs are identical
success =   all(rr1(:) == rr2(:)) && ...
            all(cc1(:) == cc2(:)) && ...
            (num_samples1 == num_samples2);
        
% Display result of test
if success, disp('Test succeeded'); [t1,t2]
else        error('Test failed');
end


