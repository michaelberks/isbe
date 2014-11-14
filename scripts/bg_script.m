%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for investigated modelling subtracted background patches
% 23rd November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Form structure of target regions
% 16 for now, look for more later
target_list(1).name = 'o04_004LML.bmp';
target_list(1).rows = [4198 4851]; %what do I mean by flipped?
target_list(1).cols = [1335 2141];
target_list(1).idx = 6;

target_list(2).name = 'o04_016LCC.bmp';
target_list(2).rows = [765 1477];
target_list(2).cols = [757 1517];
target_list(2).idx = 21;

target_list(3).name = 'o04_016LML.bmp';
target_list(3).rows = [1628 2318];
target_list(3).cols = [843 1512];
target_list(3).idx = 22;

target_list(4).name = 'o04_018RCC.bmp';
target_list(4).rows = [4856 5593];
target_list(4).cols = [4710 5329];
target_list(4).idx = 25;

target_list(5).name = 'o04_021RCC.bmp';
target_list(5).rows = [2257 3033];
target_list(5).cols = [3259 3961];
target_list(5).idx = 27;

target_list(6).name = 'o04_040LCC.bmp';
target_list(6).rows = [2140 2833];
target_list(6).cols = [1 666];
target_list(6).idx = 48; % Although there is 040RCCb

target_list(7).name = 'o04_040LML.bmp';
target_list(7).rows = [1586 2293];
target_list(7).cols = [380 1063];
target_list(7).idx = 50; %Altough there is 040RMLa

target_list(8).name = 'o04_056LML.bmp';
target_list(8).rows = [1918 2985];
target_list(8).cols = [435 1765];
target_list(8).idx = 77;

target_list(9).name = 'o04_075LCC.bmp';
target_list(9).rows = [4713 5260];
target_list(9).cols = [1492 2064];
target_list(9).idx = 101;

target_list(10).name = 'o04_078RML.bmp';
target_list(10).rows = [1407 2157];
target_list(10).cols = [2731 3453];
target_list(10).idx = 104;

% target_list(11).name = 'o04_086LCC.bmp';
% target_list(11).rows = [1750 2401];
% target_list(11).cols = [50 85];
% target_list(11).idx = 0; %086RML discarded as too near edge. Find new ROI

target_list(11).name = 'o04_104LCC.bmp';
target_list(11).rows = [3811 4510];
target_list(11).cols = [345 1014];
target_list(11).idx = 135;

target_list(12).name = 'o04_104LML.bmp';
target_list(12).rows = [1976 2720];
target_list(12).cols = [885 1587];
target_list(12).idx = 136;

target_list(13).name = 'o04_109RML.bmp';
target_list(13).rows = [2771 3529];
target_list(13).cols = [2980 3846];
target_list(13).idx = 140;

target_list(14).name = 'o04_112RCC.bmp';
target_list(14).rows = [2316 3192];
target_list(14).cols = [2740 3581];
target_list(14).idx = 146;

target_list(15).name = 'o04_125RCC.bmp';
target_list(15).rows = [1500 2150];
target_list(15).cols = [3280 3944];
target_list(15).idx = 158;

%%
% Check regions look sensible
for ii = 1:15
    i1 = imread(['C:\isbe\mammograms\new_CAD\BMP_2004\', target_list(ii).name]);
    figure; imagesc(i1); colormap(gray(256)); axis image; hold on;
    clear i1;
    plot([target_list(ii).cols(1) target_list(ii).cols(2)...
        target_list(ii).cols(2) target_list(ii).cols(1) target_list(ii).cols(1)],...
        [target_list(ii).rows(1) target_list(ii).rows(1)...
        target_list(ii).rows(2) target_list(ii).rows(2) target_list(ii).rows(1)]);
end
%% 
% Extract regions and save as .bmp patches (for use with CJR code)
for ii = 1:15
    i1 = imread(['C:\isbe\mammograms\new_CAD\BMP_2004\', target_list(ii).name]);
    patch = i1(target_list(ii).rows(1):target_list(ii).rows(2),...
        target_list(ii).cols(1):target_list(ii).cols(2));
    imwrite(patch, ['C:\isbe\dev\background\images\normal\normal_bg_', zerostr(ii, 3), '.bmp']);
    clear i1 patch;
end
        
