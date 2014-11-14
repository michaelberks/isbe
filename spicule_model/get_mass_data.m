function get_mass_data(mass_files, path, if_plot)

    if nargin < 3
        if_plot = 0;
    end

    N = length(mass_files);

    for ii = 1:N

        temp = load([path, mass_files(ii).name]);
        mass = temp.mass; clear temp;
        mass_name = ['m', mass_files(ii).name(3:end)];
        display([mass_name, '']);
        mass_ROI = mass.mass_ROI;
        subtract_ROI = mass.mass_sub_it;
        subtract_ROI(subtract_ROI < 0) = 0;
        mass_outline = mass.mass_outline;
        %n_spicules = length(spicules);
        nipple = mass.nipple;
        clear mass;

        %
        % First get mass data
        %
        [rows cols] = size(mass_ROI);
        
        %
        % Flip ROI in right mammograms
        if strcmpi(mass_name(8), 'R')
            mass_outline(:,1) = cols + 1 - mass_outline(:,1);
            mass_ROI = fliplr(mass_ROI);
            subtract_ROI = fliplr(subtract_ROI);
            nipple(1) = cols + 1 - nipple(1);
        end
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
        %c_to_o = [c_to_o(m_idx:end, :); c_to_o(1:m_idx-1, :)];

        if if_plot
            figure, imagesc(mass_ROI), colormap gray, axis image, hold on;
            plot(mass_outline(:,1), mass_outline(:,2));
            plot([mass_centroid(1) nipple(1)] , [mass_centroid(2) nipple(2)], 'g');
            plot(mass_centroid(1), mass_centroid(2), 'gx');
            plot(mean(mass_outline(:,1)), mean(mass_outline(:,2)), 'co');
            plot(nipple(1), nipple(2), 'rx');
            plot(mass_outline(1,1), mass_outline(1,2), 'yx');
        end

        mass.mass_outline = [mass_outline(:,1) - mass_centroid(1)...
            mass_outline(:,2) - mass_centroid(2)];
        
        mass.mass_ROI = mass_ROI;
        mass.subtract_ROI = subtract_ROI;
        %mass.n_spicules = n_spicules;
        mass.mass_length = m_lengths(end);
        mass.mass_area = mass_area;
        mass.mass_centroid = mass_centroid;
        mass.nipple = nipple;
        
        save([path, mass_name], 'mass');
    end
end