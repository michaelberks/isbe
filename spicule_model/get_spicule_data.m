function [mass_data spic_data] = ...
    get_spicule_data(mass_file, width, num_pts, if_plot)

if nargin < 4
    if_plot = 0;
end

temp = load(mass_file);
mass = temp.mass; clear temp;
spicules = mass.mass_spicules;
mass_ROI = mass.mass_ROI;
mass_outline = mass.mass_outline;
n_spicules = length(spicules);
nipple = mass.nipple;
clear mass;

%
% First get mass data
%
[rows cols] = size(mass_ROI);
mass_mask = roipoly(rows, cols, mass_outline(:,1), mass_outline(:,2));

[start_r start_c] = find(mass_mask, 1);
mass_outline = bwtraceboundary(mass_mask, [start_r, start_c], 'E');
mass_outline = [mass_outline(:,2) mass_outline(:,1)];
m_diff = diff(mass_outline);

mass_outline = mass_outline([true; m_diff(:,1) | m_diff(:,2)],:);

% Get mass centroid using region props
% Use nipple value to calculate starting point on mass
rp = regionprops(bwlabel(mass_mask, 4), 'Centroid', 'Area');
mass_centroid = rp.Centroid; mass_area = rp.Area; clear rp;
%mass_centroid = [mass_centroid(2) mass_centroid(1)];

m_diff = diff(mass_outline);
m_lengths = [0; cumsum(sqrt(m_diff(:,1).^2 + m_diff(:,2).^2))];
clear m_diff;

c_to_n = (nipple - mass_centroid) / sqrt(sum((nipple - mass_centroid).^2));
c_to_n = repmat(c_to_n, size(mass_outline, 1), 1);

c_to_o = mass_outline - repmat(mass_centroid, size(mass_outline,1), 1);
c_to_o = c_to_o ./ repmat(sqrt(sum(c_to_o'.^2))', 1, 2);

[md m_idx] = min(sum((c_to_o - c_to_n)'.^2));
mass_outline = [mass_outline(m_idx:end, :); mass_outline(1:m_idx-1, :)];
c_to_o = [c_to_o(m_idx:end, :); c_to_o(1:m_idx-1, :)];

if if_plot
    figure, imagesc(mass_ROI), colormap gray, axis image, hold on;
    plot(mass_outline(:,1), mass_outline(:,2));
    plot([mass_centroid(1) nipple(1)] , [mass_centroid(2) nipple(2)], 'g');
    plot(mass_centroid(1), mass_centroid(2), 'gx');
    plot(mean(mass_outline(:,1)), mean(mass_outline(:,2)), 'co');
    plot(nipple(1), nipple(2), 'rx');
    plot(mass_outline(1,1), mass_outline(1,2), 'yx');
end

spic_data = {};
for jj = 1:n_spicules
        
    pts = spicules(jj).outline;
    
    %Get spicule orientation
    s_vec = pts(end,:) - pts(1,:);
    m_vec = pts(1,:) - mass_centroid;
    spic_data(jj).s_orientation = atan2(s_vec(2), s_vec(1))...
        - atan2(m_vec(2), m_vec(1));
    %acos(dot(s_vec, m_vec) / sqrt(sum(s_vec.^2)*sum(m_vec.^2)));
    
    c_to_s = (pts(1,:) - mass_centroid) / sqrt(sum((pts(1,:) - mass_centroid).^2));
    c_to_s = repmat(c_to_s, size(mass_outline, 1), 1);
    [md s_idx] = min(sum((c_to_o - c_to_s)'.^2));
    
    %Get spicule boundary location
    spic_data(jj).s_location = m_lengths(s_idx) / m_lengths(end);
    
    %Get spicule distance
    spic_data(jj).s_distance = sum(m_vec.^2) / ...
        sum((mass_outline(s_idx,:) - mass_centroid).^2);
    
    s_diff = diff(pts);
    s_lengths = [0; cumsum(sqrt(s_diff(:,1).^2 + s_diff(:,2).^2))];
    clear s_diff;
    
    spicule = interp1(s_lengths, pts, [0:s_lengths(end)]);
    s_shape = interp1(s_lengths, pts, linspace(0, s_lengths(end), num_pts));
    clear s_lengths;
    
    [fx, fy] = gradient(spicule);
    fy = fy ./ [sqrt(sum((fy.^2)')'), sqrt(sum((fy.^2)')')]; clear fx;

    for ii = 1:length(fy(:,1)) %= number of rows in spicule
        n1_x = spicule(ii, 1) - width*fy(ii, 2);
        n1_y = spicule(ii, 2) + width*fy(ii, 1);
        n2_x = spicule(ii, 1) + width*fy(ii, 2);
        n2_y = spicule(ii, 2) - width*fy(ii, 1);

        if (n1_x >= 1 && n1_x < cols && n2_x >= 1 && n2_x < cols && n1_y >= 1 && n1_y < rows && n2_y >= 1 && n2_y < rows)
            profile(:, ii) = ...
                improfile(mass_ROI, [n1_x, n2_x], [n1_y, n2_y], 2*width+1);
        end
    end
   
    spic_data(jj).spicule = spicule;
    spic_data(jj).s_shape = ...
        s_shape - repmat(s_shape(1,:), num_pts, 1); clear s_shape;
    spic_data(jj).s_length = ...
        sqrt(sum((spicule(1,:) - spicule(end,:)).^2)); clear spicule;
    spic_data(jj).profile = profile; clear profile;
    
    if if_plot
        plot(pts(:,1), pts(:,2), 'm');
        plot([pts(1,1) pts(1,1)+s_vec(1)], [pts(1,2) pts(1,2)+s_vec(2)], 'r')
        plot([pts(1,1) pts(1,1)-m_vec(1)], [pts(1,2) pts(1,2)-m_vec(2)], 'g')
        plot(mass_outline(s_idx,1), mass_outline(s_idx,2), 'cx')
        %plot(spic_data(jj).s_shape(:,1), spic_data(jj).s_shape(:,2));
    end
    
end

mass_data.mass_outline = [mass_outline(:,1) - mass_centroid(1) mass_outline(:,2) - mass_centroid(2)];
mass_data.mass_ROI = mass_ROI;
mass_data.n_spicules = n_spicules;
mass_data.mass_length = m_lengths(end);
mass_data.mass_area = mass_area;
mass_data.mass_centroid = mass_centroid;
mass_data.nipple = nipple;