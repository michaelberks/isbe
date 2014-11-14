%Script of work for initial attempts at implementing a random forest
%builder

%Load in the zip code dataset
training_data = load('C:\isbe\dev\classification\data\zip.train');
%
profile on
tic; 
rf_mb = mb_random_forest_class_train('X', training_data(:,2:end), 'y', training_data(:,1),...
    'ntrees', 500, 'mvars', 10, 'do_oob_trace', 1);
toc;
profile off
profsave(profile('info'),'C:\isbe\dev\classification\random_forests\code_profiles\my_rf_code_zip')
whos
%
figure('Name','MB random forest OOB error rate');
plot(rf_mb.oob_trace); title('MB random forest OOB error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
save C:\isbe\dev\classification\random_forests\code_profiles\my_rf rf_mb
%clear rf_mb
%%
%Profile the Breiman/Wang implementation for building 50 trees
profile on
tic; rf_brei = classRF_train(training_data(:,2:end),training_data(:,1),500,10); toc;
profile off
profsave(profile('info'),'C:\isbe\dev\classification\random_forests\code_profiles\breiman_rf_code_zip')
whos
%
save C:\isbe\dev\classification\random_forests\code_profiles\breiman rf_brei
clear rf_brei
%%
%%

%load the twonorm dataset 
load classification\randomforest-matlab\rF_Class_C\data\twonorm    
 
%modify so that training data is NxD and labels are Nx1, where N=#of
%examples, D=# of features

X = inputs';
Y = outputs;

[N D] =size(X);
%randomly split into 250 examples for training and 50 for testing
randvector = randperm(N);

X_trn = X(randvector(1:250),:);
Y_trn = Y(randvector(1:250));
X_tst = X(randvector(251:end),:);
Y_tst = Y(randvector(251:end));
%
display('Breiman random forest: 500 trees');
tic; model = classRF_train(X_trn,Y_trn); toc;
figure('Name','Breiman random forest OOB error rate');
plot(model.errtr(:,1)); title('Breiman random forest OOB error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
Y_fit = classRF_predict(X_tst,model);
fprintf('Breiman random forest test error rate %f\n',   mean(Y_fit~=Y_tst));
%%
display('MB bagging: 500 trees');
tic; [random_forest] = mb_random_forest_class_train('X', X_trn, 'y', Y_trn, 'ntrees', 500, 'mvars', 0, 'do_oob_trace', 1); toc;
figure('Name','MB bagging OOB error rate');
plot(random_forest.oob_trace); title('MB bagging OOB error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
Y_fit = str2double(mb_random_forest_class_predict(random_forest, X_tst));
fprintf('bagging test error rate %f\n',   mean(Y_fit~=Y_tst));


display('MB random forest: 500 trees');
tic; [random_forest] = mb_random_forest_class_train('X', X_trn, 'y', Y_trn, 'ntrees', 500, 'mvars', 4, 'do_oob_trace', 1); toc;
figure('Name','MB random forest OOB error rate');
plot(random_forest.oob_trace); title('MB random forest OOB error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
Y_fit = str2double(mb_random_forest_class_predict(random_forest, X_tst));
fprintf('MB random forest test error rate %f\n',   mean(Y_fit~=Y_tst));


