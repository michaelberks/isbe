function ok = display_array(ARRAY)

% rescale the array to lie between 0 and 255;
minval = min(ARRAY);
minval = min(minval);
maxval = max(ARRAY);
maxval = max(maxval);

ARRAY=ARRAY-minval;
ARRAY=ARRAY*255/(maxval-minval);

% cast back to uint8 (0-255)
ARRAY=uint8(ARRAY);
imshow(ARRAY);
ok=1;
