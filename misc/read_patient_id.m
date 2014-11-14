function [] = read_patient_id(tag, template_dir, if_plot)
%READ_PATIENT_ID_BATCH *Insert a one line summary here*
%   [] = read_patient_id_batch()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 13-May-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if nargin < 3
    if_plot = 0;
end
if if_plot
    figure; imagesc(tag); axis image; colormap(gray(256)); hold on;
end

letters = 'ABCEGILMNOPRSTUY1234';
num_letters = length(letters);

%Get size of tag image
[rows cols] = size(tag);

%Match big letters and numbers
%%%

%Pre-allocate space for correlations cores
corr_scores = zeros([rows cols num_letters]);

%loop through each big letter
for ii = 1:num_letters
    
    %Load in big letter template
    template = u_load([template_dir letters(ii) '_big.mat']);
    [rows_t cols_t] = size(template);
    rd2 = floor(rows_t/2);
    cd2 = floor(cols_t/2);
    
    %Compute correlation scores
    c_temp = normxcorr2(template, tag);
    corr_scores(:,:,ii) = c_temp(rd2+(1:rows), cd2+(1:cols));

end

%take maximum over each template
[big_corr_max big_corr_idx] = max(corr_scores, [], 3);

%Set up mask of where we expect to find big letters numbers - this is a
%quick (hardcoded) hack for now but can be refined...
letter_mask = false(rows, cols);
letter_mask(1:100,:) = 1;
letter_mask(:,1100:end) = 1;

%Get local maxima - check more carefully size of the exclusion radius (set
%to 15 - approx 3/4 of letter width)
[big_maxima_pos] = local_image_maxima(big_corr_max, 15, letter_mask, 0.7);

%Match little numbers 
%%%
for ii = 0:9
     %Load in little number template
    template = u_load([template_dir num2str(ii) '.mat']);
    [rows_t cols_t] = size(template);
    rd2 = floor(rows_t/2);
    cd2 = floor(cols_t/2);
    
    %Compute correlation scores
    c_temp = normxcorr2(template, tag);
    corr_scores(:,:,ii+1) = c_temp(rd2+(1:rows), cd2+(1:cols));
end

%take maximum over each template
[little_corr_max little_corr_idx] = max(corr_scores, [], 3);

%Get local maxima - check more carefully size of the exclusion radius (set
%to 15 - approx 3/4 of letter width) - also which just using the inverse of
%the big letter mask, but again we could do something more sophisticated to
%specify where we expect to find little numbers)
[little_maxima_pos] = local_image_maxima(little_corr_max, 10, ~letter_mask, 0.7);

if if_plot
    %Loop through big letter maxima and plot on image of tag
    for ii = 1:size(big_maxima_pos,1)
        x = big_maxima_pos(ii,1);
        y = big_maxima_pos(ii,2);
        text(x-10, y, letters(big_corr_idx(y, x)),...
            'color', 'r', 'fontsize', 28);
    end
    
    %Loop through little maxima and plot on image of tag
    for ii = 1:size(little_maxima_pos,1)
        x = little_maxima_pos(ii,1);
        y = little_maxima_pos(ii,2);
        text(x-5, y, num2str(little_corr_idx(y, x)-1),...
            'color', 'r', 'fontsize', 16);
    end
end

%Other stuff to do:

%Work out where rows of text are from y coordinates of big/little maxima
%are - discard maxima that are outliers to these rows

%Sort maxima position by x-coordinate in each row, and read off stream of
%numbers/letters

%Decide what it is you want to return as output

%Things to consider:
% - Work out where spaces occur in Names (e.g. look for large gaps in
% x-coordinates)
% - Search for 'PMA' specifically?
% - check that date formats come out as expected
