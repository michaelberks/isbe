function [lo_roi hi_roi] = lo_v_hi_mag_ncm(lo_im, hi_im)

f1 = figure('WindowStyle', 'normal', 'units', 'normalized', 'position', [0 0 1 1]); imgray(lo_im);
f2 = figure('WindowStyle', 'normal', 'units', 'normalized', 'position', [0 0 1 1]); imgray(hi_im);

figure(f1);
[x, y, ~] = impixel();
lo_roi = lo_im(y+(-240:239), x+(-640:639), 2);
clf;
imgray(lo_roi);

figure(f2);
[x, y, ~] =impixel();
hi_roi = hi_im(y+(-240:239), x+(-640:639), 1);
clf;

close(f1);
close(f2);