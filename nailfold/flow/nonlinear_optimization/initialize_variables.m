function [flow1, flow2, displacements] = ...
    initialize_variables(gt, wt_gt, wt_noise, f_use_uv)

if (nargin==0 && nargout==0), test(); return; end

if ~exist('f_use_uv','var'), f_use_uv = true; end;

if f_use_uv
    % Zeros, plus some ground truth and some noise
    uu = zeros(size(gt.flow)) + ...
         wt_gt * real(gt.flow) + ...
         wt_noise * randn(size(gt.flow)) * 0.1;

    vv = zeros(size(gt.flow)) + ...
         wt_gt * imag(gt.flow) + ...
         wt_noise * randn(size(gt.flow)) * 0.1;
     
    flow1 = uu; flow2 = vv;
else
    % Zeros, plus some ground truth and some noise
    F = zeros(size(gt.flow)) + ...
        wt_gt * abs(gt.flow) + ...
        wt_noise * rand(size(gt.flow));
    
    theta = zeros(size(gt.flow)) + ...
            wt_gt * angle(gt.flow) + ...
            wt_noise * rand(size(gt.flow))*2*pi;

%     theta = angle(gt.orientation);
    F = ones(size(gt.orientation));
    theta = zeros(size(gt.orientation));
%     F = zeros(size(gt.orientation));

    flow1 = F; flow2 = theta;
end

displacements = zeros(size(gt.displacements)) + ...
                wt_gt * gt.displacements + ...
                wt_noise * randn(size(gt.displacements)) * 0.1;



function test()
clc;
