function [filename,window_width] = get_calibfilm_info(file_num)

filenames = zeros(50,12);
filenames = char(filenames);

filenames(1,:) = 'swcaf001.bmp';
filenames(2,:) = 'swcaf002.bmp';
filenames(3,:) = 'swcaf003.bmp';
filenames(4,:) = 'swcaf004.bmp';
filenames(5,:) = 'swcaf005.bmp';
filenames(6,:) = 'swcaf006.bmp';
filenames(7,:) = 'swcaf007.bmp';
filenames(8,:) = 'swcaf008.bmp';
filenames(9,:) = 'swcaf009.bmp';
filenames(10,:) = 'swcaf010.bmp';
filenames(11,:) = 'swcaf011.bmp';
filenames(12,:) = 'swcaf012.bmp';
filenames(13,:) = 'swcaf013.bmp';
filenames(14,:) = 'swcaf014.bmp';
filenames(15,:) = 'swcaf015.bmp';
filenames(16,:) = 'swcaf016.bmp';
filenames(17,:) = 'swcaf017.bmp';
filenames(18,:) = 'swcaf018.bmp';
filenames(19,:) = 'swcaf019.bmp';
filenames(20,:) = 'swcaf020.bmp';
filenames(21,:) = 'swcaf021.bmp';
filenames(22,:) = 'swcaf022.bmp';
filenames(23,:) = 'swcaf023.bmp';
filenames(24,:) = 'swcaf024.bmp';
filenames(25,:) = 'swcaf025.bmp';
filenames(26,:) = 'swcaf026.bmp';
filenames(27,:) = 'swcaf027.bmp';
filenames(28,:) = 'swcaf028.bmp';
filenames(29,:) = 'swcaf029.bmp';
filenames(30,:) = 'swcaf030.bmp';
filenames(31,:) = 'swcaf031.bmp';
filenames(32,:) = 'swcaf032.bmp';
filenames(33,:) = 'swcaf033.bmp';
filenames(34,:) = 'swcaf034.bmp';
filenames(35,:) = 'swcaf035.bmp';
filenames(36,:) = 'swcaf036.bmp';
filenames(37,:) = 'swcaf037.bmp';
filenames(38,:) = 'swcaf038.bmp';
filenames(39,:) = 'swcaf039.bmp';
filenames(40,:) = 'swcaf040.bmp';
filenames(41,:) = 'swcaf041.bmp';
filenames(42,:) = 'swcaf042.bmp';
filenames(43,:) = 'swcaf043.bmp';
filenames(44,:) = 'swcaf044.bmp';
filenames(45,:) = 'swcaf045.bmp';
filenames(46,:) = 'swcaf046.bmp';
filenames(47,:) = 'swcaf047.bmp';
filenames(48,:) = 'swcaf048.bmp';
filenames(49,:) = 'swcaf049.bmp';
filenames(50,:) = 'swcaf050.bmp';
filenames(51,:) = 'swcaf051.bmp';
filenames(52,:) = 'swcaf052.bmp';
filenames(53,:) = 'swcaf053.bmp';
filenames(54,:) = 'swcaf054.bmp';
filenames(55,:) = 'swcaf055.bmp';
filenames(56,:) = 'swcaf056.bmp';
filenames(57,:) = 'swcaf057.bmp';
filenames(58,:) = 'swcaf058.bmp';
filenames(59,:) = 'swcaf059.bmp';
filenames(60,:) = 'swcaf060.bmp';
filenames(61,:) = 'swcaf061.bmp';
filenames(62,:) = 'swcaf062.bmp';
filenames(63,:) = 'swcaf063.bmp';
filenames(64,:) = 'swcaf064.bmp';
filenames(65,:) = 'swcaf065.bmp';
filenames(66,:) = 'swcaf066.bmp';
filenames(67,:) = 'swcaf067.bmp';
filenames(68,:) = 'swcaf068.bmp';
filenames(69,:) = 'swcaf069.bmp';
filenames(70,:) = 'swcaf070.bmp';
filenames(71,:) = 'swcaf071.bmp';
filenames(72,:) = 'swcaf072.bmp';
filenames(73,:) = 'swcaf073.bmp';
filenames(74,:) = 'swcaf074.bmp';
filenames(75,:) = 'swcaf075.bmp';
filenames(76,:) = 'swcaf076.bmp';
filenames(77,:) = 'swcaf077.bmp';
filenames(78,:) = 'swcaf078.bmp';
filenames(79,:) = 'swcaf079.bmp';
filenames(80,:) = 'swcaf080.bmp';
filenames(81,:) = 'swcaf081.bmp';
filenames(82,:) = 'swcaf082.bmp';

window_widths = zeros(50);

window_widths(1) = 2774;
window_widths(2) = 756;
window_widths(3) = 2728;
window_widths(4) = 2448;
window_widths(5) = 2420;
window_widths(6) = 2428;
window_widths(7) = 772;
window_widths(8) = 3016;
window_widths(9) = 1908;
window_widths(10) = 2396;
window_widths(11) = 3812;
window_widths(12) = 3816;
window_widths(13) = 3812;
window_widths(14) = 3808;
window_widths(15) = 3804;
window_widths(16) = 3808;
window_widths(17) = 3792;
window_widths(18) = 3772;
window_widths(19) = 3756;
window_widths(20) = 3816;
window_widths(21) = 3816;
window_widths(22) = 3812;
window_widths(23) = 3812;
window_widths(24) = 3808;
window_widths(25) = 3804;
window_widths(26) = 3800;
window_widths(27) = 3788;
window_widths(28) = 3564;
window_widths(29) = 3780;
window_widths(30) = 3796;
window_widths(31) = 3804;
window_widths(32) = 3797;
window_widths(33) = 3797;
window_widths(34) = 3797;
window_widths(35) = 3797;
window_widths(36) = 3797;
window_widths(37) = 468;
window_widths(38) = 492;
window_widths(39) = 2760;
window_widths(40) = 1856;
window_widths(41) = 2712;
window_widths(42) = 3068;
window_widths(43) = 2264;
window_widths(44) = 3816;
window_widths(45) = 3808;
window_widths(46) = 3800;
window_widths(47) = 3788;
window_widths(48) = 3804;
window_widths(49) = 3804;
window_widths(50) = 3804;
window_widths(51) = 3816;
window_widths(52) = 3812;
window_widths(53) = 3812;
window_widths(54) = 3772;
window_widths(55) = 3776;
window_widths(56) = 3796;
window_widths(57) = 3804;
window_widths(58) = 3808;
window_widths(59) = 3812;
window_widths(60) = 3816;
window_widths(61) = 3812;
window_widths(62) = 3816;
window_widths(63) = 3800;
window_widths(64) = 3776;
window_widths(65) = 3792;
window_widths(66) = 3804;
window_widths(67) = 3808;
window_widths(68) = 3812;
window_widths(69) = 3812;
window_widths(70) = 3816;
window_widths(71) = 3816;
window_widths(72) = 3796;
window_widths(73) = 3768;
window_widths(74) = 3788;
window_widths(75) = 3796;
window_widths(76) = 3804;
window_widths(77) = 3804;
window_widths(78) = 3812;
window_widths(79) = 3812;
window_widths(80) = 3808;
window_widths(81) = 3816;
window_widths(82) = 3820;

filename = filenames(file_num,:);
window_width = window_widths(file_num);