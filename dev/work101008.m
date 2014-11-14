R = 32;
R2 = round(1.3*R);

xx = repmat(1:2*R2, 2*R2, 1) - R2;
yy = xx';
template = R^2 - xx.^2 - yy.^2;
template(xx.^2 + yy.^2 > R^2) = 0;

template2 = template;
template2(xx.^2 + yy.^2 > R2^2) = NaN;

mammo = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_024RML.mat');
mammo = imresize(mammo, 0.5, 'bilinear');
%%
xc1 = normxcorr2(template, double(mammo));
xc1 = xc1(R2:end-R2,R2:end-R2);

xc2 = mb_normxcorr2(template2, double(mammo));
figure; imagesc(template); axis image;
%
figure; imagesc(xc1); axis image; colormap(jet(256));
figure; imagesc(xc2); axis image; colormap(jet(256));
%%
m_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*mat');
mkdir C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\
for ii = 1:length(m_list)
    display(['processing mammogram: ' num2str(ii) ' of ' num2str(length(m_list))]);
    mammo = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\' m_list(ii).name]);
    original_size = size(mammo);
    mammo = imresize(mammo, 0.5, 'bilinear');

    xc_scales = -ones(size(mammo));
    scales = uint8(zeros(size(mammo)));
    for R = 5*(2:18)
        R2 = round(1.3*R);

        xx = repmat(1:2*R2, 2*R2, 1) - R2;
        yy = xx';
        template = R^2 - xx.^2 - yy.^2;
        template(xx.^2 + yy.^2 > R^2) = 0;

        xc1 = normxcorr2(template, double(mammo));
        xc1 = xc1(R2:end-R2,R2:end-R2);
        %figure; imagesc(xc1); axis image; colormap(jet(256));

        xc_scales(xc1 > xc_scales) = xc1(xc1 > xc_scales);
        scales(xc1 > xc_scales) = R;
    end
    xc_scales = imresize(xc_scales, original_size, 'bilinear');
    save(...
        ['C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\' m_list(ii).name(1:end-4) '_mass.mat'],...
        'xc_scales', 'scales');
    clear xc_scales scales mammo;
end
%%
m_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*mat');
mkdir C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\scales\
for ii = 1:length(m_list)
    load(...
        ['C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\' m_list(ii).name(1:end-4) '_mass.mat'],...
        'xc_scales', 'scales');
    load(...
        ['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\' m_list(ii).name(1:end-4) '_mask.mat']);
    mask = imresize(mask, size(xc_scales));
    
    xc_scales(~mask) = 0;
    save(...
        ['C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\' m_list(ii).name(1:end-4) '_mass.mat'],...
        'xc_scales');
    save(...
        ['C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\scales\' m_list(ii).name(1:end-4) '_mass.mat'],...
        'scales');
end
%%
m_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*mat');
mkdir C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\combined\
for ii = 1:length(m_list)
    try
        mass_map = u_load(...
            ['C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\' m_list(ii).name(1:end-4) '_mass.mat']);
        rad_map = u_load(...
            ['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' m_list(ii).name(5:end-4) '_rad_map_064.mat']);

        combined_map = rad_map .* mass_map;

        save(...
            ['C:\isbe\asymmetry_project\data\mass_maps\2004_screening\abnormals\combined\' m_list(ii).name(1:end-4) '_mass.mat'],...
            'combined_map');
    catch
        display(['Bugger, ' m_list(ii).name(1:end-4) ' didn''t work']);
    end
end
    