%Load in the zip code dataset
training_data = load('C:\isbe\dev\classification\data\zip.train');
%%
profile on
tic; 
rf_mb = mb_random_forest_class_train('X', training_data(y,2:end), 'y', training_data(y,1),...
    'ntrees', 20, 'mvars', 10, 'do_oob', 0);
toc;
profile off
profsave(profile('info'),'C:\isbe\dev\classification\random_forests\code_profiles\my_rf_code_zip')
whos
%%
figure('Name','MB random forest OOB error rate');
plot(rf_mb.oob_trace); title('MB random forest OOB error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
save C:\isbe\dev\classification\random_forests\code_profiles\my_rf rf_mb
%%
y = randsample(1:size(training_data,1),100000, true)'; 
%Profile the Breiman/Wang implementation for building 50 trees
profile on
tic; rf_brei = classRF_train(training_data(y,2:end),training_data(y,1),20,10); toc;
profile off
profsave(profile('info'),'C:\isbe\dev\classification\random_forests\code_profiles\breiman_rf_code_zip')
whos
%
save C:\isbe\dev\classification\random_forests\code_profiles\breiman rf_brei
clear rf_brei
%%
x = sortrows([X(noderows,jvar) Cnode], 1); % get sorted jth x variable
            
             % Determine if there's anything to split along this variable
             maxeps = max(eps(x(1,1)), eps(x(end,1)));
             if x(1,1)+maxeps > x(end,1)
                continue;
             end
             rows = find(x(1:end-1,1)+maxeps < x(2:end,1));
             if isempty(rows)
                continue;
             end

            %Ccum = cumsum(Cnode(idx,:));         % cum. class counts
            Ccum = cumsum(x(:,2));
%%
tic; [x1 idx] = sortrows(training_data, 2); toc;

[r c] = size(training_data);
tic;
%x2 = zeros(r,c);
[irank, ifail] = m01de(training_data, int32(1), int32(2), 'A', 'n2', int32(2));
m01ea(x(:,ii), int32(1), irank);
% for ii = 1:1
%     x2(:,c) = m01ea(x(:,ii), int32(1), irank);
% end
%%
display('sort and rank matlab')
tic
[x1 idx] = sort(training_data(:,2));
toc
%
display('sort nag')
tic
[irank, ifail] = m01da(training_data(:,2), int32(1), 'A');
%x2 = m01ea(training_data(:,2), int32(1), irank);
x2 = training_data(irank,2);
toc;
%
display('sort nag')
tic
x3 = m01ca(training_data(:,2), int32(1), 'A');
toc;
%%
profile on
tic; rf_brei_bar = classRF_train(training_data, training_labels,100,4); toc;
profile off
profsave(profile('info'),'C:\isbe\dev\classification\random_forests\code_profiles\breiman_rf_code_bar')
whos
%
save C:\isbe\dev\classification\random_forests\code_profiles\breiman rf_brei_bar
clear rf_brei_bar
%%
bar = create_gauss_bar(10, 10, 45, 64, 64);


figure;
subplot(1,2,1); mesh(bar);
subplot(1,2,2); mesh(bar_interp1);
%%
t1 = zeros(10,1);
t2 = zeros(10,1);
for ii = 1:10

    tic; bar_interp1 = interp2(1:64, (1:64)', bar, 1:.5:64, (1:.5:64)', 'linear'); t1(ii) = toc;
    tic;
    [px, py, lamda, mu, c, ifail] = e01da(1:64, 1:64, bar(:));
    bar_interp2 = reshape(e02df(1:.5:64, (1:.5:64)', lamda, mu, c), [], length(1:.5:64));
    t2(ii) = toc;
end
%%
t1 = zeros(10,1);
t2 = zeros(10,1);
for goat = 1:10
    tic;
    [px, py, lamda_r, mu_r, c_r, ifail] = e01da(x(1,:), y(:,1), real(z_unwrap(:)));
    [px, py, lamda_i, mu_i, c_i, ifail] = e01da(x(1,:), y(:,1), imag(z_unwrap(:)));
    z_out_unwrap = reshape(complex(e02df(xi(1,:), yi(:,1), lamda_r, mu_r, c_r), e02df(xi(1,:), yi(:,1), lamda_i, mu_i, c_i)), size(xi));
    t2(goat) = toc;
    
    tic; z_out_unwrap = interp2(x, y, z_unwrap, xi, yi, 'bicubic'); t1(goat) = toc;
end


%%
figure;
subplot(1,2,1); mesh(bar);
subplot(1,2,2); mesh(bar_interp2);
%%
mb_dual_tree_transform