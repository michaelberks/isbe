root_dir = 'F:\Disney Movies\';
mkdir([root_dir '/segments']);
movie_list = dir([root_dir '*.avi']);

for i_m = 2:length(movie_list)
    input_avi = ['"' root_dir movie_list(i_m).name '"'];
    for i_s = 1:12
        output_avi = ['"' root_dir 'segments/' movie_list(i_m).name(1:end-4) zerostr(i_s,2) '.avi"'];
        st = zerostr((i_s-1)*5,2);
        cmd = ...
            ['ffmpeg -i ' input_avi ' -vcodec copy -acodec copy -ss 00:' st ':00 -t 00:05:00 ' output_avi];
        [status,result] = system(cmd) %#ok
        %display(cmd);
    end
end
        
        