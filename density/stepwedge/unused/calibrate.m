iptsetpref('ImshowBorder','tight'); % show less border round image
iptsetpref('TruesizeWarning','off'); % don't clutter up output with image size warnings

%clear
%disp('using correction image already in memory');

IMCORR = imread('IMCORR_LARGE_2004.BMP','bmp');

% request calibration image name
%file_num = input('what calibration film number? ');
    fid = fopen('calibresults.txt','aw');

for file_num=1:1,
    
    [filenames,window_width] = get_calibfilm_info_2004(file_num);
    
    disp(['reading in film: ' filenames]);
    
    % read in image to be calibrated
    [IMAGE,map] = imread(filenames);
    IMAGE = imresize(IMAGE,0.2);
    
    if(1)
        
        % CIMAGE = curve_corr(IMAGE,map,IMCORR,2,0,window_width);
        CIMAGE = geomcorr(IMAGE,IMCORR,1,window_width);
        
    end
    
    %CIMAGE = IMAGE;
    
    % next need to call wedgecal
    
    disp('First select small stepwedge area:');
    wedgevals1 = wedgecal_new(CIMAGE,map,1);
    disp('Now select large stepwedge area:');
    wedgevals2 = wedgecal_new(CIMAGE,map,2);
    
    % plot the wedgevals to prove it worked:
    figure(4);
    plot(wedgevals1);
    hold on
    plot(wedgevals2,'r');
    hold off
    
    % create lookup table to convert pixel value to x_sw
    sw_lookup1 = sw_lookup_table(wedgevals1);
    sw_lookup2 = sw_lookup_table(wedgevals2);
    figure(5);
    plot(sw_lookup1);
    hold on
    plot(sw_lookup2,'r');
    hold off
    
    disp('Now select ROI in centre of phantom: ');
    figure(6);
    mask = roipoly(CIMAGE);
    
    phantom_val = sum(sum(double(mask).*double(CIMAGE)))/sum(sum(double(mask)))
    
    small_step = sw_lookup1(round(phantom_val))
    large_step = sw_lookup2(round(phantom_val))
    
    fprintf(fid,'%i %7.4f %7.4f %7.4f \n',file_num,phantom_val,small_step,large_step);
end
    fclose(fid);
