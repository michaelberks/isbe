function [m_ii covar_ii sample_sorted] =...
    eda_gaussian(fn, no_samples, no_it, m_init, covar_init, sample_ii)
%
% The skeleton of this function performs a pdf estimation based
% optimisation algorithm - however it has been corrupted to use
% specifically with optimising mass_model weights and no longer works for
% generalised functions.

plot_no = no_it - 20;

no_samples = no_samples + rem(no_samples, 2); %ensures even sample size
sample_half = no_samples / 2; 
dim = length(m_init);

if nargin < 6
    sample_ii = sample_from_normal(m_init, covar_init, no_samples);
end

temp = load('C:\isbe\dev\mass_model.mat');
mass_model = temp.mass_model; clear temp;
    
for ii = 1:no_it,
    f = zeros(no_samples, 1);
    for jj = 1:no_samples
        xy = sample_ii(jj, :);
        f(jj) = mean(feval(fn, mass_model, xy));
    end
    sample_f = [sample_ii, f]; clear sample_ii;
    sample_sorted = sortrows(sample_f, dim+1); %sort the sample in ascending order of f values
    
    if ii > plot_no,
        if dim == 2
            figure('WindowStyle', 'docked')
            hold on;
            plot(sample_sorted(1:sample_half, 1),...
                sample_sorted(1:sample_half, 2), 'g.');
            plot(sample_sorted(sample_half + 1:end, 1),...
                sample_sorted(sample_half + 1:end, 2), 'r.');
            xlabel('shape'); ylabel('texture');
        elseif dim == 3
            figure('WindowStyle', 'docked')
            hold on;
            plot3(sample_sorted(1:sample_half, 1), ...
                sample_sorted(1:sample_half, 2),...
                sample_sorted(1:sample_half, 3), 'g.');
            plot3(sample_sorted(sample_half + 1:end, 1),...
                sample_sorted(sample_half + 1:end, 2),...
                sample_sorted(sample_half + 1:end, 3), 'r.');
            xlabel('shape'); ylabel('texture'); zlabel('scale');
        end
        
    end
    
    m_ii = mean(sample_sorted(1:sample_half, 1:dim));
    covar_ii = cov(sample_sorted(1:sample_half, 1:dim));
    sample_ii = sample_from_normal(m_ii', covar_ii, no_samples);
end
    
