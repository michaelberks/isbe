function cost = match_cost(ref, tgt, method, f_translate, f_scale)
if (nargin==0 && nargout==0), test_script(); return; end

if ~exist('method','var'), method = 'normxcorr2'; end
if ~exist('f_translate','var'), f_translate = false; end
if ~exist('f_scale','var'), f_scale = false; end

% Transform the reference patch if requested
if (f_scale)
    ref = ref - mean(ref(:));
    sd = sqrt(ref(:)'*ref(:) / (numel(ref)-1));
    ref = ref / sd;
elseif (f_translate)
    ref = ref - mean(ref(:));
end

switch method
    case {'normxcorr2'},
        cost = normxcorr2(double(ref), double(tgt));
        patch_hw = (size(ref)-1) / 2;
        cost = cost((1+patch_hw(1)):(end-patch_hw(1)), ...
                    (1+patch_hw(2)):(end-patch_hw(2)));
        % minimize negative cost to maximize correlation
        cost = -cost; 
                
    case {'filter2', 'xcorr2'}
        cost = filter2(double(ref), double(tgt));
        cost = cost(end:-1:1, end:-1:1);
        % minimize negative cost to maximize correlation
        cost = -cost; 
                
    case {'phasecorr'}
        % Phase correlation goes wonky when the reference patch is
        % normalized to zero mean: the first component of the FFT,
        % corresponding to the DC part (i.e. the mean), is zero such that
        % the first element of the product and its absolute value are also
        % zero, giving a divide-by-zero error.
        
        border = (size(tgt)-size(ref)) / 2;
        ref2 = mb_pad(ref, border, 0);
        prd = fft2(double(ref2)) .* conj(fft2(double(tgt)));
        cost = fftshift(ifft2(prd ./ abs(prd)));
        % minimize negative cost to maximize correlation
        cost = -cost; 

    case {'sse'},
        if (f_scale)
            cost = sse_scale(ref, tgt);
        elseif (f_translate)
            cost = sse_translate(ref, tgt);
        else
            cost = sse_no_transform(ref, tgt);
        end
        
    case {'mse'}
        cost = func(ref, tgt, 'sse', f_translate, f_scale);
        cost = cost / numel(ref);

    case {'sabs'},
        patch_sz = size(ref);
        patch_hw = (patch_sz - 1) / 2;
        tgt2 = mb_pad(tgt, patch_hw, NaN);
        ni = patch_sz(1)-1;
        nj = patch_sz(2)-1;
        cost = nan(size(tgt));
        for i = 1:size(tgt,1)
            for j = 1:size(tgt,2)
                subtgt = tgt2(i:(i+ni), j:(j+nj));
                if (f_scale)
                    subtgt = subtgt - mean(subtgt(:));
                    sd = sqrt(subtgt(:)'*subtgt(:) / (numel(subtgt)-1));
                    subtgt = subtgt / sd;
                elseif (f_translate)
                    subtgt = subtgt - mean(subtgt(:));
                end
                cost(i,j) = sum(abs(ref(:) - subtgt(:)));
            end
        end

    case {'mabs'}
        cost = func(ref, tgt, 'sabs', f_translate, f_scale);
        cost = cost / numel(ref);
end

function cost = sse_scale(ref, tgt, patch_sz)
patch_sz = size(ref);
patch_hw = (patch_sz - 1) / 2;
tgt2 = mb_pad(tgt, patch_hw, NaN);
cost = nan(size(tgt));
for i = 1:size(tgt,1)
    for j = 1:size(tgt,2)
        subtgt = tgt2(i:(i+patch_sz(1)-1), j:(j+patch_sz(2)-1));
        
        subtgt = subtgt - mean(subtgt(:));
        sd = sqrt(subtgt(:)'*subtgt(:) / (numel(subtgt)-1));
        subtgt = subtgt / sd;
        
        diff = ref - subtgt;
        cost(i,j) = diff(:)'*diff(:);
    end
end

function cost = sse_translate(ref, tgt, patch_sz)
patch_sz = size(ref);
patch_hw = (patch_sz - 1) / 2;
tgt2 = mb_pad(tgt, patch_hw, NaN);
cost = nan(size(tgt));
for i = 1:size(tgt,1)
    for j = 1:size(tgt,2)
        subtgt = tgt2(i:(i+patch_sz(1)-1), j:(j+patch_sz(2)-1));
        
        subtgt = subtgt - mean(subtgt(:));

        diff = ref - subtgt;
        cost(i,j) = diff(:)'*diff(:);
    end
end

function cost = sse_no_transform(ref, tgt, patch_sz)
patch_sz = size(ref);
patch_hw = (patch_sz - 1) / 2;
tgt2 = mb_pad(tgt, patch_hw, NaN);
ni = patch_sz(1)-1;
nj = patch_sz(2)-1;
cost = nan(size(tgt));
for i = 1:size(tgt,1)
    for j = 1:size(tgt,2)
        diff = ref - tgt2(i:(i+ni), j:(j+nj));
        cost(i,j) = diff(:)' * diff(:);
    end
end


%% Test script
function test_script()
clc;

c = 5;
x = randn(1, (2*c)+1);
w = 1;
c = c + 1; % Offset the centre point
t = x((c+1-w):(c+1+w));

mc = func(t, x, 'normxcorr2')
mc = func(t, x, 'xcorr2')
mc = func(t, x, 'filter2')
mc = func(t, x, 'phasecorr')
mc = func(t, x, 'mse')

f_profile = true;

% Check timing
if (f_profile)
    profile off;
    profile clear;
    profile on;
end
tic;
for i = 1:20
    ref = rand(21,21);
    tgt = rand(101,101);
    cost = func(ref, tgt, 'sse');
%     cost = func(ref, tgt, 'normxcorr2');
end
toc;
if (f_profile)
    profile off;
    profile report;
end

