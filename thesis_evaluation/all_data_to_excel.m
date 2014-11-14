function all_data_to_excel(num_users, excel_file)

all_ratings = NaN(90, num_users);
all_timings = NaN(90, num_users);
all_ids = cell(1,num_users);
for ii = 1:num_users
    user_id = ['User' zerostr(ii,2)];
    load(['C:\isbe\dev\observer_study\user_data\user', zerostr(ii,2)]);
    user_data_to_excel(user_data, user_id, excel_file);
    for jj = 1:60;
        all_ratings(user_data.ratings(jj,1),ii) = user_data.ratings(jj,2);
        all_timings(user_data.ratings(jj,1),ii) = user_data.ratings(jj,3);
    end
    all_ids{ii} = user_id;    
end
xlswrite(excel_file, {'Mass Index'}, 'All ratings', 'A1');
xlswrite(excel_file, all_ids, 'All ratings', 'B1');
xlswrite(excel_file, (1:90)', 'All ratings', 'A2');
xlswrite(excel_file, all_ratings, 'All ratings', 'B2');

xlswrite(excel_file, {'Mass Index'}, 'All timings', 'A1');
xlswrite(excel_file, all_ids, 'All timings', 'B1');
xlswrite(excel_file, (1:90)', 'All timings', 'A2');
xlswrite(excel_file, all_timings, 'All timings', 'B2');
