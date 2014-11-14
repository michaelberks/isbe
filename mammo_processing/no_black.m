function no_black(input_file, output_file, safe)

if nargin < 3
    safe = 1;
end
    
if strcmp(input_file, output_file) && safe
    display('error: cant overwrite file when safe is on.'); 
else    
    image = imread(input_file);

    [x, y] = size(image);
    for i=1:x,
        for j=1:y,
            if image(i, j) == 0,
                image(i, j) = 1;
            end
        end
    end
    imwrite(image, output_file);
    clear image;
end