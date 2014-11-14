% SAMPLE_NEW_SPICULES1 - sample new spicules from spicule appearance model 
%
% Usage: [new_spicules] = sample_new_spicules1(spicule_AM, n_spicules,...
%                                               width, if_plot)
% 
% Arguments:
% spicule_AM - spicule appearance model
% n_spicules - scalar integer, number of new spicules to sample
% width - scalar integer, the half width of the new spicule normal profiles
% if_plot - [0 1] turn on/off displayinf new spciules
%
% Returns: 
% new_spicules - n_spicules-length structure of new spicules
%
% Notes:
% 
%
% See also: GENERATE_SPICULE_AM SAMPLE_NEW_MASSES
%
% References:
%
% Author:   Michael Berks
%           Imaging Science and Biomedical Engineering
%           University of Manchester
%
function [new_spicules] = sample_new_spicules1(spicule_AM, n_spicules, width, if_plot);

    mean_s  = spicule_AM.mean_s;
    P_s     = spicule_AM.P_s;
    B_s     = spicule_AM.B_s;
    L_s     = spicule_AM.L_s;
    mean_w  = spicule_AM.mean_w;
    P_w     = spicule_AM.P_w;
    B_w     = spicule_AM.B_w;
    L_w     = spicule_AM.L_w;
    mean_b  = spicule_AM.mean_b;
    P_b     = spicule_AM.P_b;
    B_b     = spicule_AM.B_b;
    L_b     = spicule_AM.L_b;
    mean_p  = spicule_AM.mean_p;
    P_p     = spicule_AM.P_p;
    B_p     = spicule_AM.B_p;
    L_p     = spicule_AM.L_p;
    mean_l  = spicule_AM.mean_l;
    P_l     = spicule_AM.P_l;
    B_l     = spicule_AM.B_l;
    L_l     = spicule_AM.L_l;

    mean_c  = spicule_AM.mean_c;
    P_c     = spicule_AM.P_c;
    B_c     = spicule_AM.B_c;
    L_c     = spicule_AM.L_c;

    W_s     = spicule_AM.W_s;
    W_w     = spicule_AM.W_w;
    W_b     = spicule_AM.W_b;
    W_p     = spicule_AM.W_p;
    W_l     = spicule_AM.W_l;

    m = length(mean_s) / 2;

    k_s = length(L_s);
    k_w = length(L_w);
    k_b = length(L_b);
    k_p = length(L_p);
    
    %
    % Generate specified number of masses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ii = 1;
    while ii <= n_spicules

        % sample new combined appearance vectors - assume normal distribution
        % of modes
        %%%

        % Compute new shape vector, texture vector and scale
        new_c = (randn(length(L_c), 1) .* sqrt(L_c));

        Q_s = P_c(1:k_s,:); 
        Q_w = P_c(k_s+1:k_s + k_w,:);
        Q_b = P_c(k_s+k_w+1:k_s+k_w+k_b,:);
        Q_p = P_c(k_s+k_w+k_b+1:end-1,:);
        Q_l = P_c(end, :);

        new_s = mean_s + (P_s*Q_s*new_c)' / W_s;
        new_w = mean_w + (P_w*Q_w*new_c)' / W_w;
        new_b = mean_b + (P_b*Q_b*new_c)' / W_b;
        new_p = mean_p + (P_p*Q_p*new_c)' / W_p;
        new_l = mean_l + (P_l*Q_l*new_c)' / W_l;
        
        new_b(new_b < 0) = 0;
        new_w(new_w < 1) = 1;

        new_shape(:,1) = new_s(1:m);
        new_shape(:,2) = new_s(m+1:2*m);
        new_shape = new_shape*new_l / sqrt(sum((new_shape(end,:) - new_shape(1,:)).^2));
        
        clear new_s;
        new_spicules(ii).new_shape = new_shape;
        %new_spicules(ii).s_orientation = new_p(1);
        %new_spicules(ii).s_location = new_p(2);
        %new_spicules(ii).s_distance = new_p(3);
        new_spicules(ii).s_orientation = 0.25*randn(1);
        new_spicules(ii).s_location = rand(1);
        new_spicules(ii).s_distance = 1 + 0.1*randn(1);
        
        %calculate length of new spicule
        V = diff(new_shape);
        N = [0; cumsum(sqrt(V(:,1).^2 + V(:,2).^2))];
        
        %create landmark points spaced along pixel length for new shape and
        %profile
        n_pts = round(N(end) / 5);%
        
        if n_pts <= 1, continue; end;
        
        % Reconstruct profile appearance vectors
        profile_x = [1:N(end)+1];
        prows = 2*width + 1; pcols = length(profile_x);
        
        %calculate widths at landmarks along profile from new_w
        w_vector = interp1(linspace(1, pcols, length(new_w)), new_w, 1:pcols);
        w_matrix = repmat(w_vector, prows, 1);
        
        %calculate brightness along profile from new_b
        b_vector = interp1(linspace(1, pcols, length(new_b)), new_b, 1:pcols);
        b_matrix = repmat(b_vector, prows, 1)*0.25;

        offsets = repmat([-width:width]', 1, pcols);
        offsets(abs(offsets) >= w_matrix) = w_matrix(abs(offsets) >= w_matrix);
        
        new_spicules(ii).profile =...
            b_matrix .* sqrt(1 - (offsets.^2 ./ w_matrix.^2));
        
        new_spicules(ii).landmarks = ...
            interp1(N, new_shape, linspace(0, N(end), n_pts));
          
        new_spicules(ii).landmarks_w = interp1(linspace(1, pcols,...
            length(new_w)), new_w, linspace(1, pcols, n_pts));
        
        ii = ii + 1;
    end
end