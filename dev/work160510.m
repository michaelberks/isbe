load C:\isbe\misc\mnist\mnist_train.mat
load C:\isbe\misc\mnist\mnist_test.mat
%%
mkdir C:\isbe\misc\mnist\images;
mkdir C:\isbe\misc\mnist\images\train;
mkdir C:\isbe\misc\mnist\images\test;

%%
f_tr = fopen('C:\isbe\misc\mnist\images\train_labels.txt', 'wt');
for i_tr = 1:size(train_X,1)
    imwrite(reshape(train_X(i_tr,:), 28, 28)',...
        ['C:\isbe\misc\mnist\images\train\im' zerostr(i_tr,5) '.png']);
    
    fprintf(f_tr, '%s %d\n', ['im' zerostr(i_tr,5) '.png'], train_labels(i_tr));
end
fclose(f_tr);

f_te = fopen('C:\isbe\misc\mnist\images\test_labels.txt', 'wt');
for i_te = 1:size(test_X,1)
    imwrite(reshape(test_X(i_te,:), 28, 28)',...
        ['C:\isbe\misc\mnist\images\test\im' zerostr(i_te,5) '.png']);
    
    fprintf(f_te, '%s %d\n', ['im' zerostr(i_te,5) '.png'], test_labels(i_te));
end
fclose(f_te);
%%
num_ims = 60000;
train_X = zeros(num_ims, 28*28, 'uint8');
train_X2 = zeros(num_ims, 28*28, 'uint8');

for i_im = 1:num_ims
    im1 = imread(['C:\isbe\misc\mnist\images\train\im' zerostr(i_im,5) '.png']);
    im2 = imread(['C:\isbe\misc\mnist\images\train2\im' zerostr(i_im-1,5) '.png']);
    
    train_X(i_im,:) = im1(:)';
    train_X2(i_im,:) = im2(:)';
end
%%
train_X = double(train_X);
train_X2 = double(train_X2);

available = true(num_ims,1);
all_idx = 1:num_ims;
idx2 = zeros(num_ims,1);
for i_im = 1:num_ims
    local_idx = find(~any(bsxfun(@minus, train_X2(available,:), train_X(i_im,:)),2));
    if isempty(local_idx)
        display(['Warning, no match found for im ' num2str(i_im)]);
    elseif length(local_idx) > 1
        display(['Warning, im ' num2str(i_im) ' matches multiple images: ' num2str(local_idx(:)')]);
    else
        available_idx = all_idx(available);
        global_idx = available_idx(local_idx);
        idx2(i_im) = global_idx;
        available(global_idx) = 0;
    end
end
    