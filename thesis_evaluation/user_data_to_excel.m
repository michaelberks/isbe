function user_data_to_excel(user_data, user_id, excel_file)

xlswrite(excel_file, {'Mass Index'}, user_id, 'A1');
xlswrite(excel_file, {'Mass Rating'}, user_id, 'B1');
xlswrite(excel_file, {'Time (s)'}, user_id, 'C1');
xlswrite(excel_file, {'Feedback'}, user_id, 'D1');

xlswrite(excel_file, user_data.ratings, user_id, 'A2');
if isfield(user_data, 'feedback');
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    n = size(user_data.feedback, 2)-1;
    syn_idx = user_data.ratings(:,2) <= 2;
    feedback = NaN(60, n);
    
    feedback(syn_idx,:) = user_data.feedback(:,2:end);
    feedback_other = cell(60,1);
    feedback_other(syn_idx) = user_data.feedback_other;
    xlswrite(excel_file, feedback, user_id, 'D2');
    xlswrite(excel_file, feedback_other, user_id, [alphabet(4+n), '2']);
end
