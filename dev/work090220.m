% Script for transferring IPl phase information from mass region into a
% normal region

%Long term, we need to choose which masses are really suitable for this -
%and ideally, only choose the area containing the abnormalities associated
%with the mass

%First we need to resize our mass regions so they are divisible by 6
load('C:\isbe\dev\files\bg_files.mat');
mkdir C:\isbe\dev\background\images\mass_2_16
%
for ii = 1:length(bg_files);
    
    load(['C:\isbe\dev\masses\', bg_files(ii).name]);
    mass_bg = double(mass.mass_ROI) - mass.subtract_ROI;
    
    %Reduce resolution by a half
    mass_bg = imresize(mass_bg, 0.5, 'bilinear');
    [r c] = size(mass_bg);
    if r < 512
        rem_r = rem(r, 16);
        r_start = ceil(rem_r / 2) + 1;
        r_end = r - floor(rem_r / 2);
    else
        r_start = ceil((r - 512) / 2) + 1;
        r_end = r - floor((r - 512) / 2);
    end
    if c < 512
        rem_c = rem(c, 16);
        c_start = ceil(rem_c / 2) + 1;
        c_end = c - floor(rem_c / 2);
    else
        c_start = ceil((c - 512) / 2) + 1;
        c_end = c - floor((c - 512) / 2);
    end
        
    mass_bg = mass_bg(r_start:r_end, c_start:c_end);
    
    imwrite(uint8(mass_bg), ['C:\isbe\dev\background\images\mass_2_16\mass_bg_', zerostr(ii,3), '.bmp']);
    
end
%%
% Now build dual_trees of the new images
mb_build_dual_tree(...
    'ImageDir', 'C:\isbe\dev\background\images\mass_2_16\',...
    'OutputDir', 'C:\isbe\dev\background\dual_tree\mass_2_16\');

%%
% Now try putting the mass_bg ILP + ICP coefficients in a normal background
% - it'll be very interesting to see the effect of masss with non-breast
% areas in (e.g. mass001)
mass_list = dir('C:\isbe\dev\background\dual_tree\mass_2_16\*.mat');
normal_list = dir('C:\isbe\dev\background\dual_tree\normal512\*.mat');

idx = round(length(normal_list)*rand);
normal_tree = u_load(['C:\isbe\dev\background\dual_tree\normal512\', normal_list(idx).name]);
[normal_ilp normal_icp] = mb_dual_tree_transform(normal_tree);
normal_region = dtwaveifm2(normal_tree);
%%
for ii = [3 4 6 8 16 24 27 34 36 40 51 52 53]
    mass_tree = u_load(['C:\isbe\dev\background\dual_tree\mass_2_16\', mass_list(ii).name]);
    [mass_ilp mass_icp] = mb_dual_tree_transform(mass_tree);
    
    [r c] = size(mass_tree{4}(:,:,1));
    rs = ceil((32 - r) / 2);
    cs = ceil((32 - c) / 2);
    re = 32 - floor((32 - r) / 2);
    ce = 32 - floor((32 - c) / 2);
    
    
    combined_icp = normal_icp(1:3);
    combined_ilp = normal_ilp(1:3);
    
    combined_tree1 = normal_tree;
    for lev = 1:4
        
        rs_lev = rs * 2^(4-lev);
        cs_lev = cs * 2^(4-lev);
        re_lev = re * 2^(4-lev);
        ce_lev = ce * 2^(4-lev);
        
        if lev < 4
            combined_icp{lev}(rs_lev+1:re_lev, cs_lev+1:ce_lev,:) = mass_icp{lev};
            combined_ilp{lev}(rs_lev+1:re_lev, cs_lev+1:ce_lev,:) = mass_ilp{lev};
        end
        combined_tree1{lev}(rs_lev+1:re_lev, cs_lev+1:ce_lev,:) = mass_tree{lev};
    end
    combined_icp{4} = [];
    combined_ilp{4} = combined_tree1{4};
    combined_tree = mb_dual_tree_transform_i(combined_ilp, combined_icp);   
    combined_tree(5:6) = normal_tree(5:6);%#ok 
    
    for lev = 1:4
        combined_tree{lev} = abs(combined_tree1{lev}) .* exp(i*angle(combined_tree{lev}));
    end
    
    combined_region = dtwaveifm2(combined_tree);
    combined_region1 = dtwaveifm2(combined_tree1);
    mass_region = dtwaveifm2(mass_tree);
    figure; 
%     subplot(1,3,1); imagesc(normal_region); axis image; colormap(gray(256));
%     subplot(1,3,2); imagesc(mass_region); axis image; colormap(gray(256));
%     subplot(1,3,3); imagesc(combined_region1); axis image; colormap(gray(256));
imagesc(combined_region1); axis image; colormap(gray(256));
%     figure; 
%     subplot(1,2,1); imagesc(combined_region); axis image; colormap(gray(256));
%     subplot(1,2,2); imagesc(combined_region1); axis image; colormap(gray(256));
end

%%
% I'm concerned - why are do the mass regions have lower local variance
% than the normal regions. Is this because they way they have been
% downsampled is different?

%Don't worry, this is sorted now, see generating_regions_script. Next...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%First let look at the 4th scaling coefficients to see if swapping between
%them is feasible - e.g. it might not be