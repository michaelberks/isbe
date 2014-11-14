%**************************************************************************
%***** Script to produce results for CVPR submission Nov 2011 using *******
%****************** the DRIVE retinography dataset ************************
%**************************************************************************
%**************************************************************************
warning('off', 'load_uint8:missing_variables');

clear;

%% 1. Produce response, orientation and scale maps for the analytic methods
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\training\predictions\orientation\'];
do_mono =0;
do_g1d = 0;
do_g2d = 1;

if do_mono
    mkdir([retroot 'mono\analytic\orientations']);
    mkdir([retroot 'mono\analytic\responses']);
    mkdir([retroot 'mono\analytic\scales']);
end
if do_g1d
    mkdir([retroot 'g1d\analytic\orientations']);
    mkdir([retroot 'g1d\analytic\responses']);
    mkdir([retroot 'g1d\analytic\scales']);
end
if do_g2d
    mkdir([retroot 'g2d\analytic\orientations']);
    mkdir([retroot 'g2d\analytic\orientations_c']);
    mkdir([retroot 'g2d\analytic\responses']);
    mkdir([retroot 'g2d\analytic\scales']);
end

for ii = 21:40
    %load retinogram and merge RGB channels
    ret = u_load([asymmetryroot('shared'),'data\retinograms\DRIVE\training\images_extended\' zerostr(ii,2) '_training_ext.mat']);
    ret = rgb2gray(ret);
    
    if do_mono
        %Compute repsonses for mono and save
        [response_map d ori_map scale_map] = monogenic_multiscale(ret, 4, 4, 2, 0.65);
        save_uint8([retroot 'mono\analytic\orientations\' zerostr(ii,3) '_training_ori.mat'], ori_map);
        save_uint8([retroot 'mono\analytic\responses\' zerostr(ii,3) '_training_response.mat'], response_map);
        save_uint8([retroot 'mono\analytic\scales\' zerostr(ii,3) '_training_scale.mat'], scale_map);
    end
    if do_g1d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_1st_derivative_gradient2(ret, [1 2 4 8]);
        save_uint8([retroot 'g1d\analytic\orientations\' zerostr(ii,3) '_training_ori.mat'], ori_map);
        save_uint8([retroot 'g1d\analytic\responses\' zerostr(ii,3) '_training_response.mat'], response_map);
        save_uint8([retroot 'g1d\analytic\scales\' zerostr(ii,3) '_training_scale.mat'], scale_map);
    end
    if do_g2d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_clover_line(ret, [1 2 4 8]);
        save_uint8([retroot 'g2d\analytic\orientations\' zerostr(ii,3) '_training_ori.mat'], ori_map);
        save_uint8([retroot 'g2d\analytic\orientations_c\' zerostr(ii,3) '_training_ori.mat'], complex(cos(2*ori_map), sin(2*ori_map)));
        save_uint8([retroot 'g2d\analytic\responses\' zerostr(ii,3) '_training_response.mat'], response_map);
        save_uint8([retroot 'g2d\analytic\scales\' zerostr(ii,3) '_training_scale.mat'], scale_map);
    end
    

end