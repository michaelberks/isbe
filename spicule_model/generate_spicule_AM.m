% GENERATE_SPICULE_AM - Top level function to generate appearance model
%
% Usage: [spicule_model] = generate_spicule_AM(spicule_data)
% 
% Arguments:
% spicule_data - N-length structure, where N is the number of spicules to
% be modeled, containing 5 fields:
%  w_vector - vector of parameters defining the spicule width
%  b_vector - vector of parameters defining the spicule brightness
%  o_vector - vector of offsets from the centre along the spicule length           
%  s_vector - m*2 matrix of (x y)-co-ordinates, where m is the number of 
%   landmark points, definig the spicule shape 
%  s_length - scalar of spicule length 
% This structure is created automatically using get_spicule_parameters
%
% Returns:      spicule_model - 
%  
% Notes:
% 
%
% See also: GENERATE_MASS_AM GET_SPICULE_PARAMETERS
%
% References:
%
% Author:   Michael Berks
%           Imaging Science and Biomedical Engineering
%           University of Manchester
%
function [spicule_AM, w_matrix, s_matrix] = generate_spicule_AM(spicule_data)
    
    colors = ['r' 'g' 'b' 'y' 'm' 'k' 'c'];
    N = length(spicule_data); %N = number of spicules in the model
    m = size(spicule_data(1).s_vector, 1);
    %
    % Build shape model of spicules
    %

    % Put spicule shape, width, brightness vectors into matrices
    
    %f1 = figure; f2 = figure;
    for ii = 1:N
        sv = spicule_data(ii).s_vector;
        th = -atan2(sv(end, 2), sv(end, 1));
        
        %figure(f1); plot(sv(:,1), sv(:,2), colors(mod(ii,7)+1)); hold on;
        sv = [cos(th) -sin(th); sin(th) cos(th)] * sv';
        %figure(f2); plot(sv(1,:), sv(2,:), colors(mod(ii,7)+1)); hold on;
        
        s_matrix(ii,:) = [sv(1,:), sv(2,:)];
        w_matrix(ii,:) = spicule_data(ii).w_vector;
        b_matrix(ii,:) = spicule_data(ii).b_vector;
        %o_matrix(ii,:) = spicule_data(ii).o_vector;
        s_lengths(ii) = sqrt(sum((sv(:,1) - sv(:,end)).^2));%spicule_data(ii).s_length;
        p_matrix(ii,1) = spicule_data(ii).s_orientation;
        p_matrix(ii,2) = spicule_data(ii).s_location;
        p_matrix(ii,3) = spicule_data(ii).s_distance;
    end
    
    [mean_l, P_l, B_l, L_l] = pca(s_lengths', 0.98); %compute principal modes of scale

    sf = repmat((mean_l./s_lengths)', 1, 2*m); %make matrix of scale factors
    s_matrix = s_matrix .* sf; %scale shape matrix
    
%     figure; hold on;
%     for ii = 1:N
%         
%         plot(s_matrix(ii, 1:end/2), s_matrix(ii, (end/2) + 1:end), colors(mod(ii,7)+1));
%     end
    
    test_widths = polyval(w_matrix(1,:), s_matrix(1, 1:end/2));
    [s_matrix(1, 1:end/2)' test_widths'] 
    [mean_s, P_s, B_s, L_s] = pca(s_matrix, 0.98); %compute principal modes of shape
    [mean_w, P_w, B_w, L_w] = pca(w_matrix, 0.98); %compute principal modes of width
    [mean_b, P_b, B_b, L_b] = pca(b_matrix, 0.98); %compute principal modes of brightness
    [mean_p, P_p, B_p, L_p] = pca(p_matrix, 0.98); %compute principal modes of position
    
    %
    %Calculate weights for combined model
    %%%%%%%%%%%%%%%%%%%%%%
    k_s = length(L_s);
    k_w = length(L_w);
    k_b = length(L_b);
    k_p = length(L_p);
    
    W_s = k_s / sum(sqrt(L_s)); %sqrt()
    W_w = k_w / sum(sqrt(L_w));
    W_b = k_b / sum(sqrt(L_b));
    W_p = k_p / sum(sqrt(L_p));
    W_l = 1 / sqrt(L_l); %length(L_l = 1)

    combined_data = [W_s*B_s; W_w*B_w; W_b*B_b; W_p*B_p; W_l*B_l]';

    [mean_c, P_c, B_c, L_c] = pca(combined_data, 0.98);
    
    spicule_AM.mean_s  = mean_s;
    spicule_AM.P_s     = P_s;
    spicule_AM.B_s     = B_s;
    spicule_AM.L_s     = L_s;
    spicule_AM.mean_w  = mean_w;
    spicule_AM.P_w     = P_w;
    spicule_AM.B_w     = B_w;
    spicule_AM.L_w     = L_w;
    spicule_AM.mean_b  = mean_b;
    spicule_AM.P_b     = P_b;
    spicule_AM.B_b     = B_b;
    spicule_AM.L_b     = L_b;
    spicule_AM.mean_p  = mean_p;
    spicule_AM.P_p     = P_p;
    spicule_AM.B_p     = B_p;
    spicule_AM.L_p     = L_p;
    spicule_AM.mean_l  = mean_l;
    spicule_AM.P_l     = P_l;
    spicule_AM.B_l     = B_l;
    spicule_AM.L_l     = L_l;

    spicule_AM.mean_c  = mean_c;
    spicule_AM.P_c     = P_c;
    spicule_AM.B_c     = B_c;
    spicule_AM.L_c     = L_c;

    spicule_AM.W_s     = W_s;
    spicule_AM.W_w     = W_w;
    spicule_AM.W_b     = W_b;
    spicule_AM.W_p     = W_p;
    spicule_AM.W_l     = W_l;
    
end
