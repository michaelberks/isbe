
function CIMAGE = geomcorr(IMAGE,IMCORR,leftright,window_width)

% flip correction image if this is a left breast
if(leftright)
    disp('flipping IMCORR');
    IMCORR=flipdim(IMCORR,1);
    IMCORR=flipdim(IMCORR,2);
end

% resize IMCORR to be same size as IMAGE if necessary
dimIMAGE = size(IMAGE);
dimIMCORR = size(IMCORR);
if(dimIMAGE(1)==dimIMCORR(1) & dimIMAGE(2)==dimIMCORR(2))
    IMCORR_USE = IMCORR;
else
    IMCORR_USE = imresize(IMCORR,dimIMAGE);
end

IMCORR_USE = uint8(double(IMCORR_USE)*3797/window_width);
CIMAGE = uint8(max(0,min(255,double(IMAGE) + double(IMCORR_USE))));

