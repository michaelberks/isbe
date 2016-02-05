load('C:\isbe\nailfold\data\2_year_study\results\miccai\auto_stats.mat');
load('C:\isbe\nailfold\data\2_year_study\results\miccai\people_stats.mat');

%Check what images we have each subject at each timepoint
present_mask = squeeze(sum(sum(people_stats.present,2),3));

%As long as no values in the present mask are greater than one, we can use
%the same trick of summing across 2nd and 3rd dimensions to get 2D matrix
%for each measurement
any(present_mask(:) > 1)

gradeable_mask = squeeze(sum(sum(people_stats.gradeable,2),3)); 
distal1_mask = squeeze(sum(sum(people_stats.num_distal_vessels>0,2),3));
distal2_mask = squeeze(sum(sum(people_stats.num_distal_vessels>1,2),3));

%Valid subjects are ones with images at v1, v5 and v6 (0, 1 and 2 years)
present_0_12_24 = all(present_mask(:,[1 5 6]),2);
present_0_24 = all(present_mask(:,[1 6]),2);

gradeable_0_24 = all(gradeable_mask(:,[1 6]),2);
distal1_0_24 = all(distal1_mask(:,[1 6]),2);
distal2_0_24 = all(distal2_mask(:,[1 6]),2);



%Recategorise
recats = people_stats.category;
recats(people_stats.category == 6) = 1;
recats(people_stats.category == 4) = 2;
recats(people_stats.category == 2) = 4;
recats(people_stats.category == 3) = 3;
%recats(people_stats.category == 5) = 5;

%Display numbers per group
cat_names = {'HC', 'PR', 'lSSc', 'dSSc',  'UCTD'};
cat_idx = false(length(people_stats.category),5);
cat_numbers = zeros(5,5);
for i_cat = 1:5
    cat_idx_i = recats == i_cat;
    cat_numbers(1,i_cat) = sum(cat_idx_i);
    cat_numbers(2,i_cat) = sum(cat_idx_i & present_0_12_24);
    cat_numbers(3,i_cat) = sum(cat_idx_i & present_0_24);
    cat_numbers(4,i_cat) = sum(cat_idx_i & distal1_0_24);
    cat_numbers(5,i_cat) = sum(cat_idx_i & distal2_0_24);
    display([cat_names{i_cat} ': ' num2str(cat_numbers(:,i_cat)') ]);
    
    cat_idx(:,i_cat) = cat_idx_i;
end
%, 'dSSc', 'lSSc', 
cat_colors = 'rgbym';
%%
temp_cats = auto_stats.category;
temp_cats(auto_stats.category == 6) = 1;
temp_cats(auto_stats.category == 4) = 2;
temp_cats(auto_stats.category == 2) = 4;
temp_cats(auto_stats.category == 3) = 3;
xlswrite('C:\isbe\nailfold\2_year_cats.xlsx', temp_cats, 1, 'A1');