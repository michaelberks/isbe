%%
ii = 99;
go_on = 1;
while go_on;
    ii = ii + 1;
    file_out = ['C:\isbe\dev\masses\left_out_models\model', zerostr(ii, 3)];
    one_out_files = m_files;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
    
    c = clock;
    %go_on = ~(c(3)==18 && c(4) == 9)&& ii < 50;
    go_on = ii < 100;
end
%%