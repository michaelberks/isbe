fid = fopen('C:/isbe/nailfold/data/aam/candidates/template_matching/aligned/model_qualities.txt');
q_txt = textscan(fid, '%s %f', 'delimiter', ':');
fclose(fid);
model_q = q_txt{2}; clear q_txt;
[sorted_model_qualities qidx] = sort(model_q);

feature = 'orig';

num_rows = 3;
num_cols = 4;
ii = 1;

while ii <= 338
    figure;
    for row = 1:num_rows, 
        for col = 1:num_cols
            if ii > 338
                break;
            end
            jj = qidx(ii);
            
            load(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\apex' zerostr(jj, 4) '.mat'], 'apex_candidate');
            f1 = fopen(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\out_points\orig\candidate_apex' zerostr(jj,4) '.pts']);
            textscan(f1, '%[^{]');
            fgetl(f1);
            vessel_str = textscan(f1,'%[^}]');
            fclose(f1);
            test_pts = str2num(vessel_str{1}{1});

            im = imread(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\images\orig\candidate_apex' zerostr(jj,4) '.png']);

            axes('units', 'normalized', 'position', [(col-1)/num_cols (row-1)/num_rows 1/num_cols 1/num_rows]);
            
            [r c] = size(im);
            roi = im(round(r/4):round(3*r/4), round(c/4):round(3*c/4));
            
            imagesc(im); axis image off; colormap(gray(256)); caxis([min(roi(:)) max(roi(:))]);
            hold on;
            if isempty(apex_candidate.true_vessel_xy)
                plot(test_pts(:,1), test_pts(:,2), 'bx');
            else
                plot(apex_candidate.true_vessel_xy(:,1),apex_candidate.true_vessel_xy(:,2), 'rx'); 
                plot(test_pts(:,1), test_pts(:,2), 'gx');
            end
            
            ii = ii + 1;
            
        end
    end
    
end
%%
candidate_labels = false(338,1);
for ii = 1:338
    load(['C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\apex' zerostr(ii, 4) '.mat'], 'apex_candidate');
    candidate_labels(ii) = ~isempty(apex_candidate.true_vessel_xy);
end
save('C:\isbe\nailfold\data\aam\candidates\template_matching\aligned\candidate_labels.mat', 'candidate_labels');
%%