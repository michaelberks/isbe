%%
%
% Analysing the errors from varying lengths of shape and texture vector
%

% Load errors and put means of error type in matrices
er_size_total = zeros(10, 6);
er_size_shape = zeros(10, 6);
er_size_tex = zeros(10, 6);
er_size_scale = zeros(10, 6);
er_size_combined = zeros(10, 6);

for ii = 1:10
    for jj = 1:6
        iii = 100*ii; jjj = 5000*jj;
        load(['C:\isbe\dev\size\er', num2str(iii), '_', num2str(jjj)])
        
        er_size_total(ii,jj) = mean(er.ind_error.total);
        er_size_shape(ii,jj) = mean(er.ind_error.shape);
        er_size_tex(ii,jj) = mean(er.ind_error.tex);
        er_size_scale(ii,jj) = mean(er.com_error.scale);
        er_size_combined(ii,jj) = mean(er.com_error.total);
        clear er;
    end
end

save C:\isbe\dev\size\er_size er_size*
%%
%
% Plot line graphs of the results
%
figure; plot(er_size_total);
title('Total individual model error (|M - M|) for varying lengths of shape and texture vector');
xlabel('Length of shape vector');
ylabel('Error');
legend({'5000', '10000', '15000', '20000', '25000', '30000'});
%%
figure; plot(er_size_shape);
title('Shape model error for varying lengths of shape and texture vector');
xlabel('Length of shape vector');
ylabel('Error');
legend({'5000', '10000', '15000', '20000', '25000', '30000'});

figure; plot(er_size_tex);
title('Texture model error for varying lengths of shape and texture vector');
xlabel('Length of shape vector');
ylabel('Error');
legend({'5000', '10000', '15000', '20000', '25000', '30000'});

figure; plot(er_size_scale);
title('Scale model error for varying lengths of shape and texture vector');
xlabel('Length of shape vector');
ylabel('Error');
legend({'5000', '10000', '15000', '20000', '25000', '30000'});
        