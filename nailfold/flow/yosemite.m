% PROBABLY JUNK NOW

yosroot = 'U:\data\optical_flow\yosemite';

nr = 252;
nc = 316;
nf = 14;

imgStack = zeros([nr, nc, nf]);
flowStack = zeros([nr, nc, nf]);
for i = 1:nf
    % Read image
    filename = sprintf('yos-images/yos%i/data', i-1);
    fid = fopen(fullfile(yosroot, filename), 'r');
    if fid ~= 0
        imgdata = fread(fid, inf, 'uchar');
        imgStack(:,:,i) = reshape(imgdata, [nc, nr])';
        fclose(fid);
    end

    % Read flow
    filename = sprintf('yos-flows/actual-flow%i/descriptor', i-1);
    fid = fopen(fullfile(yosroot, filename), 'r');
    if fid ~= 0
        s = fscanf(fid, '%s');
        sind = findstr(s, '(scale') + 6;
        [scalestr, tok] = strtok(s(sind:end), ')');
        scale = str2double(scalestr);
        [pedestalstr, tok] = strtok(tok(11:end), ')');
        pedestal = str2double(pedestalstr);
        fclose(fid);
    end
    filename = sprintf('yos-flows/actual-flow%i/data0', i-1);
    fid = fopen(fullfile(yosroot, filename), 'r');
    if fid ~= 0
        vdata = (double(fread(fid, inf, 'uchar')) * scale) + pedestal;
        vdata = reshape(vdata, [nc, nr])';
        fclose(fid);
    end
    filename = sprintf('yos-flows/actual-flow%i/data1', i-1);
    fid = fopen(fullfile(yosroot, filename), 'r');
    if fid ~= 0
        udata = (double(fread(fid, inf, 'uchar')) * scale) + pedestal;
        udata = reshape(udata, [nc, nr])';
        fclose(fid);
    end
    flowStack(:,:,i) = complex(udata, vdata);
end

figure(1); clf; hold off; colormap(gray(256));
for i = 1:nf
    subplot(2,1,1);
        image(uint8(imgStack(:,:,i)));
        axis('image');
    subplot(2,1,2);
        flo = cat(3, real(flowStack(:,:,i)), imag(flowStack(:,:,i)));
        image(flowToColor(flo));
        axis('image');
    pause(0.1);
end

return

% Define pixel positions in image 2
[xx2,yy2] = meshgrid(1:nc, 1:nr);

% Derive corresponding positions in image 1
xx1 = xx2 - real(flowStack(:,:,1));
yy1 = yy2 - imag(flowStack(:,:,1));

imginterp2 = interp2(imgStack(:,:,1), xx1, yy1);

figure(2); clf; hold off; colormap(gray(256));
subplot(1,2,1); 
    image(uint8(imgStack(:,:,2)));
    axis('image');
subplot(1,2,2); 
    image(uint8(imginterp2));
    axis('image');
    
figure(3); clf; hold off; colormap(jet(256));
    imgdiff = imgStack(:,:,2) - imginterp2;
    imagesc(abs(imgdiff));
    axis('image');
    disp(sum(abs( imgdiff(~isnan(imgdiff)) )));

