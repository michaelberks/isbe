function [filenames,window_width] = get_calibfilm_info(file_num)

filenames = zeros(20,14);
filenames = char(filenames);
window_widths = zeros(20);

% these are 2004 calibration films, digitised in July 2004
filename(1,:) = 'S0430G40F1.BMP';
filename(2,:) = 'S0435G15F1.BMP';
filename(3,:) = 'S0425G25F1.BMP';
filename(4,:) = 'S0420G10F1.BMP';
filename(5,:) = 'S0415G55F1.BMP';
filename(6,:) = 'S0415G35F1.BMP';
filename(7,:) = 'S0415G15F1.BMP';
filename(8,:) = 'S0410G20F1.BMP';
filename(9,:) = 'S0405G65F1.BMP';
filename(10,:) = 'S0435G15F2.BMP';
filename(11,:) = 'S0430G40F2.BMP';
filename(12,:) = 'S0430G00F2.BMP';
filename(13,:) = 'S0425G25F2.BMP';
filename(14,:) = 'S0420G10F2.BMP';
filename(15,:) = 'S0415G55F2.BMP';
filename(16,:) = 'S0415G35F2.BMP';
filename(17,:) = 'S0415G15F2.BMP';
filename(18,:) = 'S0410G20F2.BMP';
filename(19,:) = 'S0405G65F2.BMP';
filename(20,:) = 'S0405G25F2.BMP';

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
