load C:\isbe\dev\files\r_files50.mat
load C:\isbe\dev\files\r_files51.mat

for ii = 100:100:1000
    for jj = 5000:5000:50000
        file_out = ['C:\isbe\dev\models\size\model', num2str(ii), '_', num2str(jj)];
        mass_model = generate_mass_AM(...
            r_files51, file_out, ii, 'C:\isbe\dev\masses\', jj);
    end
    
end
% need to do errors, but just make models for now
%%
for ii = 200:100:1000
    for jj = 5000:5000:30000
        model_in = ['C:\isbe\dev\models\size\model', num2str(ii), '_', num2str(jj)];
        [er.com_error er.ind_error] = model_errors(model_in, r_files50, 'C:\isbe\dev\masses\');
        save(['C:\isbe\dev\size\er', num2str(ii), '_', num2str(jj)], 'er');
        clear er;
        pack;
    end   
end