function [filenames,window_width] = get_calibfilm_info(file_num)

filenames = zeros(20,14);
filenames = char(filenames);
window_widths = zeros(20);

% these are 2004 calibration films, digitised in July 2004
filenames(1,:) = 'S0430G40F1.BMP';
filenames(2,:) = 'S0435G15F1.BMP';
filenames(3,:) = 'S0425G25F1.BMP';
filenames(4,:) = 'S0420G10F1.BMP';
filenames(5,:) = 'S0415G55F1.BMP';
filenames(6,:) = 'S0415G35F1.BMP';
filenames(7,:) = 'S0415G15F1.BMP';
filenames(8,:) = 'S0410G20F1.BMP';
filenames(9,:) = 'S0405G65F1.BMP';
filenames(10,:) = 'S0435G15F2.BMP';
filenames(11,:) = 'S0430G40F2.BMP';
filenames(12,:) = 'S0430G00F2.BMP';
filenames(13,:) = 'S0425G25F2.BMP';
filenames(14,:) = 'S0420G10F2.BMP';
filenames(15,:) = 'S0415G55F2.BMP';
filenames(16,:) = 'S0415G35F2.BMP';
filenames(17,:) = 'S0415G15F2.BMP';
filenames(18,:) = 'S0410G20F2.BMP';
filenames(19,:) = 'S0405G65F2.BMP';
filenames(20,:) = 'S0405G25F2.BMP';

window_widths(1) = 4080;
window_widths(2) = 4080;
window_widths(3) = 4076;
window_widths(4) = 4080;
window_widths(5) = 4080;
window_widths(6) = 4080;
window_widths(7) = 4068;
window_widths(8) = 4052;
window_widths(9) = 4080;
window_widths(10) = 4080;
window_widths(11) = 4080;
window_widths(12) = 4080;
window_widths(13) = 4080;
window_widths(14) = 4080;
window_widths(15) = 4080;
window_widths(16) = 4080;
window_widths(17) = 4068;
window_widths(18) = 4064;
window_widths(19) = 4076;
window_widths(20) = 4040;


filenames = filenames(file_num,:);
window_width = window_widths(file_num);
