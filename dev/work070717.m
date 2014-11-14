for ii = 1:185
    
    if strcmp(files(ii).name(end-5:end-4), 'RR')
        new_str = files(ii).name;
        new_str(end-5:end-4) = [];
        movefile(files(ii).name, new_str);
    end
end