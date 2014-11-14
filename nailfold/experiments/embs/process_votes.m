test_dir = 'C:\isbe\nailfold\data\training\';
im_dir = [test_dir 'images\'];
votes_dir = [test_dir 'predictions\apex_votes\'];
apex_dir = [test_dir 'apexes\'];

sigma = 1.5;
g_width = 6;

x_o = repmat(-g_width:g_width, 2*g_width+1, 1);
y_o = x_o';
g_o = exp(-(x_o.^2 + y_o.^2)/sigma^2);
g_norm = sum(g_o(:));

num_gt = 0;
all_candidate_labels_rf = [];

for i_im = 6% 1:5
    
    apex_im = imread([im_dir 'nailfold' zerostr(i_im, 3) '.bmp']);
    [rows cols] = size(apex_im);
    vote_im = zeros(rows, cols);
    
    rows_lim = rows-g_width;
    cols_lim = cols-g_width;
    %
    votes_list = dir([votes_dir 'nailfold' zerostr(i_im, 3) '_*.txt']);
    for i_f = 1:length(votes_list)
        
        votes = load([votes_dir votes_list(i_f).name]);
        tic;
        for i_v = 1:size(votes,1)

            v_i = votes(i_v,:);

            v_ii = floor(v_i);
            v_id = v_i - v_ii;

            if (v_ii(1) < g_width) || (v_ii(2) < g_width) || (v_ii(1) > cols_lim) || (v_ii(2) > rows_lim)
                continue;
            end

            x_v = x_o - v_id(1);
            y_v = y_o - v_id(2);

            r_i1 = v_ii(2) - g_width;
            r_i2 = v_ii(2) + g_width;
            c_i1 = v_ii(1) - g_width;
            c_i2 = v_ii(1) + g_width;

            g_i = exp(-(x_v.^2 + y_v.^2)/sigma^2);

            vote_im(r_i1:r_i2,c_i1:c_i2) = vote_im(r_i1:r_i2,c_i1:c_i2) + g_i;
        end
        toc;
        clear votes;
    end
    vote_im = vote_im / g_norm;
    
    %vote_im = imfilter(vote_im, fspecial('gaussian', [10 10], 2));
    
    save([votes_dir 'nailfold' zerostr(i_im, 3) '_votes.mat'], 'vote_im');

    figure; 
    subplot(2,1,1); 
    imgray(apex_im);
    subplot(2,1,2); 
    imgray(vote_im);

end
%%




