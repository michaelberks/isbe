function arg_out = convert2cell(arg_in)
    if arg_in(1) == '{'
        str_delims = find(arg_in == '''');
        n_delims = length(str_delims);
        arg_out = cell(n_delims/2, 1);
        for ii = 2:2:n_delims
            i1 = str_delims(ii-1) + 1;
            i2 = str_delims(ii) - 1;
            arg_out{ii/2} = arg_in(i1:i2);
        end
    else
        arg_out = arg_in;
    end
end