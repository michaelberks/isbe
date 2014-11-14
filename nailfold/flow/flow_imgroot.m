function imgpath = flow_imgroot()
% Define the path where images are stored (and where output will go)

%% Use simulated images
% imgroot = 'U:\projects\nailfold\synthesis\';
%     imgpath = fullfile(imgroot, '20130729T113817'); % Slow flow
%     imgpath = fullfile(imgroot, '20130729T143620'); % Fast flow
% Idealized
%     imgpath = fullfile(imgroot, '20130902T133413'); % Steady
%     imgpath = fullfile(imgroot, '20130902T142901'); % Shaky
% Realistic
%     imgpath = fullfile(imgroot, '20130906T135325'); % Steady
%     imgpath = fullfile(imgroot, '20130902T144351'); % Shaky
% Analysed by Mike
%     imgpath = fullfile(imgroot, 'analysed/registered');
% Small set for comparison with C++ implementation
%     imgpath = fullfile(imgroot, '20130819T140643');
% d = dir(fullfile(imgroot, '*T*'));
% imgpath = fullfile(imgroot, d(end).name);

%% Use real images
imgroot = 'U:\projects\nailfold\capture\';
% 
imgpath = fullfile(imgroot, '2012_10_22\Left.Digit4.x300\');
    imgpath = fullfile(imgpath, 'seq1\corrected\registered_g1d\cropped1');
%     imgpath = fullfile(imgpath, 'seq2\corrected\registered_g1d\cropped1');
%     imgpath = fullfile(imgpath, 'seq1\corrected\registered_g1d\cropped2');
%     imgpath = fullfile(imgpath, 'seq2\corrected\registered_g1d\cropped2');

% imgpath = fullfile(imgroot, 'SRFT/short');
%     imgpath = fullfile(imgpath, 'f059/Left/Digit4/x300/cropped1');
% 
%     imgpath = fullfile(imgpath, 'f72\Left\Digit4\x300\corrected\registered_g1d\cropped1');
%     imgpath = fullfile(imgpath, 'f72\Left\Digit4\x300\corrected\registered_g1d\cropped2');
%     imgpath = fullfile(imgpath, 'f72\Left\Digit4\x300\corrected\registered_g1d\cropped3');
% 
%     imgpath = fullfile(imgpath, 'f090_v2\Left\Digit4\x300\cropped1');
%     imgpath = fullfile(imgpath, 'f090_v2\Left\Digit4\x300\cropped2');
%     imgpath = fullfile(imgpath, 'f090_v2\Left\Digit4\x300\cropped3');
%     imgpath = fullfile(imgpath, 'f090_v2\Left\Digit4\x300\cropped4');
