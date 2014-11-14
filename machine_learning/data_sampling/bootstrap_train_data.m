function [X y boot_idx] = bootstrap_train_data(X_in, y_in)

    N = size(X_in,1);
    boot_idx = ceil(N*rand(N,1));

    X = X_in(boot_idx,:);
    if nargin > 1
        y = y_in(boot_idx,:);
    end
end