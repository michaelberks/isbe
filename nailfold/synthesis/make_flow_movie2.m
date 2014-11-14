clear;

imgroot = 'U:\projects\nailfold\synthesis\showcase\ideal_example';
inpath = fullfile(imgroot, 'input');
outpath = fullfile(imgroot, 'flow');

load(fullfile(inpath,'uv.mat'));

[m,n,nFrames] = size(u_mat);

bgval = 1;

f_export = true;

u_sum = zeros(m,n);
v_sum = zeros(m,n);
nValid = zeros(m,n);
for i = 1:nFrames
    ui = u_mat(:,:,i);
    vi = v_mat(:,:,i);
    
    valid = ~isnan(ui);
    
    nValid = nValid + double(valid);
    u_sum(valid) = u_sum(valid) + ui(valid);
    v_sum(valid) = v_sum(valid) + vi(valid);
    
    rgb = complex2rgb(complex(u_sum./nValid, v_sum./nValid), ...
                      [], 1.0, [], bgval);
    rgb(isnan(rgb)) = bgval;
    
    if f_export
        filename = sprintf('frame_%04d.png', i);
        imwrite(rgb, fullfile(outpath, filename));
    else
        figure(1); clf;
            image(rgb);
            axis('image','ij');
        drawnow;
    end
end