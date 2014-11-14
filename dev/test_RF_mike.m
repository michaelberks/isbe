
A1 =[ 0.7568    0.2154    0.2032
    0.9899    0.2841    0.1344
    0.7686    0.7754    0.8108
    0.9028    0.9874    0.4549
    0.9533    0.7413    0.5824
    0.0459    0.8053    0.0997
    0.6341    0.8037    0.2238
    0.9419    0.2282    0.0819
    0.5379    0.3639    0.1084
    0.7901    0.1256    0.3808];
A2 =[ 2.1963    1.5902    2.1601
    1.5897    2.1594    1.4536
    1.4325    1.3893    1.6488
    1.3360    1.7726    1.7286
    1.2704    2.1261    1.3161
    1.4901    2.0604    2.1480
    1.2540    1.8068    1.5640
    2.1077    1.2420    1.2115
    1.9218    1.8189    1.7196
    1.3915    1.6649    2.0125];
A3 =[3.3358    3.1371    3.2564
    2.7047    3.4766    3.4943
    2.8800    3.1293    2.6008
    2.9817    2.9299    2.8854
    2.9871    2.8824    3.4194
    3.3559    3.1145    3.0396
    3.4694    2.6103    3.0930
    2.6974    2.5127    3.0148
    3.2992    2.8774    3.1355
    3.2521    3.1330    3.1384];
X_trn = [A1; A2; A3];
Y_trn = [1*ones(10,1); 2*ones(10,1); 2*ones(10,1)];

% Random Forests Breiman's code
[N nvars] = size(X_trn);
ntree=100;
random_m = ceil(sqrt(nvars));
    clear extra_options
        extra_options.importance = 1; %(0 = (Default) Don't, 1=calculate)
   extra_options.proximity = 1; %(0 = (Default) Don't, 1=calculate)
   extra_options.oob_prox = 0;
   extra_options.classwt = [1/3 2/3];
   extra_options.DEBUG_ON = true;

model = classRF_train(X_trn,Y_trn, ntree, random_m, extra_options);
    Y_hat = classRF_predict(X_trn,model);
    fprintf('\n error rate %f\n',   length(find(Y_hat~=Y_trn))/length(Y_trn));
   
% Random Forests mike's code
random_forest = mb_random_forest_class_train('X', X_trn, 'y', Y_trn, 'do_proximities', 1, 'oob_proximities', 1, 'prior', [.5 .5]);
[y_fit votes proximities] = mb_random_forest_class_predict(random_forest, X_trn);
pred_y=str2double(y_fit);
    fprintf('\n error rate %f\n',   length(find(pred_y~=Y_trn))/length(Y_trn));

