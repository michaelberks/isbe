delete('orca_nailfold_vid.avi');
aviobj = avifile('orca_nailfold_vid.avi');
aviobj.compression = 'cinepak';
aviobj.quality = 100;
aviobj.colormap = gray(236);
aviobj.fps = 9;

for ii = 1:225
    vid_frame = double(imread('N:\nailfold Capillaroscopy\oRCA Flash 2.8 Test\ORCA Flash 2.8 22ms Exp.tif', 'index', ii));
    small_frame = uint8(255*(vid_frame - min(vid_frame(:))) /(max(vid_frame(:) - min(vid_frame(:)))));
    aviobj = addframe(aviobj, small_frame);
end
aviobj = close(aviobj);
%%
vid_info = imfinfo('N:\Nailfold Capillaroscopy\ORCA-R2 Test\A Binning 2X2, normal scan.tif');
%delete('orcaR2_2_nailfold_vid.avi');
aviobj = avifile('orcaR2_2_nailfold_vid.avi');
aviobj.compression = 'cinepak';
aviobj.quality = 100;
aviobj.colormap = gray(236);
aviobj.fps = 9;

for ii = 1:length(vid_info);
    vid_frame = double(imread('N:\Nailfold Capillaroscopy\ORCA-R2 Test\A Binning 2X2, normal scan.tif', 'index', ii));
    small_frame = uint8(255*(vid_frame - min(vid_frame(:))) /(max(vid_frame(:) - min(vid_frame(:)))));
    aviobj = addframe(aviobj, small_frame);
end
aviobj = close(aviobj);
%
vid_info = imfinfo('N:\Nailfold Capillaroscopy\ORCA-R2 Test\A Binning 2X2, fast scan.tif');
%delete('orcaR2_2f_nailfold_vid.avi');
aviobj = avifile('orcaR2_2f_nailfold_vid.avi');
aviobj.compression = 'cinepak';
aviobj.quality = 100;
aviobj.colormap = gray(236);
aviobj.fps = 9;

for ii = 1:length(vid_info);
    vid_frame = double(imread('N:\Nailfold Capillaroscopy\ORCA-R2 Test\A Binning 2X2, fast scan.tif', 'index', ii));
    small_frame = uint8(255*(vid_frame - min(vid_frame(:))) /(max(vid_frame(:) - min(vid_frame(:)))));
    aviobj = addframe(aviobj, small_frame);
end
aviobj = close(aviobj);
%%

