frame_list = dir('\\isbe-matrix\transfer\For Mike\2012_08_06\mike_dmk_camera\*.png');
mkdir('\\isbe-matrix\transfer\For Mike\2012_08_06\mike_dmk_camera\enhanced\');

for ii = 0:2
    cmin = inf;
    cmax = 0;
    for jj = 1:60
        idx = 60*ii + jj;
        frame = imread(['\\isbe-matrix\transfer\For Mike\2012_08_06\mike_dmk_camera\' frame_list(idx).name]);
        cmin = min(min(frame(:)), cmin);
        cmax = max(max(frame(:)), cmax);
    end
    cmin = double(cmin);
    cmax = double(cmax);
    for jj = 1:60
        idx = 60*ii + jj;
        frame = double(imread(['\\isbe-matrix\transfer\For Mike\2012_08_06\mike_dmk_camera\' frame_list(idx).name]));
        frame = (frame - cmin) / (cmax - cmin);
        imwrite(frame, ['\\isbe-matrix\transfer\For Mike\2012_08_06\mike_dmk_camera\enhanced\' frame_list(idx).name]);
    end
end
%%

frame_list = dir('\\isbe-matrix\transfer\For Mike\2012_08_06\rohma_dmk_camera\*.png');
mkdir('\\isbe-matrix\transfer\For Mike\2012_08_06\rohma_dmk_camera\enhanced\');

for ii = 0:1
    cmin = inf;
    cmax = 0;
    for jj = 1:60
        idx = 60*ii + jj;
        frame = imread(['\\isbe-matrix\transfer\For Mike\2012_08_06\rohma_dmk_camera\' frame_list(idx).name]);
        cmin = min(min(frame(:)), cmin);
        cmax = max(max(frame(:)), cmax);
    end
    cmin = double(cmin);
    cmax = double(cmax);
    for jj = 1:60
        idx = 60*ii + jj;
        frame = double(imread(['\\isbe-matrix\transfer\For Mike\2012_08_06\rohma_dmk_camera\' frame_list(idx).name]));
        frame = (frame - cmin) / (cmax - cmin);
        imwrite(frame, ['\\isbe-matrix\transfer\For Mike\2012_08_06\rohma_dmk_camera\enhanced\' frame_list(idx).name]);
    end
end