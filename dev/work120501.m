ori_all = zeros(584, 565);
for ii = 1:20
    ori = load_uint8(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\orientation\dtg2\rf_3\' zerostr(ii,2) '_test_ext_class.mat']);
    if ismember(ii, [1 3 5 9 11 12 15]); 
        ori = sqrt(fliplr(ori));
        ori = complex(-real(ori), imag(ori)).^2;
        
    end
    figure; imgray(complex2rgb(ori));
    ori_all = ori_all + ori;
end
ori_all = ori_all / 20;
figure; imgray(complex2rgb(ori_all));
    