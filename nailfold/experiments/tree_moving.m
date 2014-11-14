tt = 1;
for ii = 1:10
    for jj = 0:9
        
        orig_name = ['C:\isbe\nailfold\models\vessel\apex_location\test_rf' num2str(ii) '_' zerostr(jj, 3) '.bfs'];
        orig_name_u = ['C:\isbe\nailfold\models\vessel\apex_location\test_rf' num2str(ii) '_' zerostr(jj, 3) '_used.bfs'];
        
        new_name = ['C:\isbe\nailfold\models\vessel\apex_location\test_rf_' zerostr(tt, 4) '.bfs'];
        new_name_u = ['C:\isbe\nailfold\models\vessel\apex_location\test_rf_' zerostr(tt, 4) '_used.bfs'];
        movefile(orig_name, new_name);
        movefile(orig_name_u, new_name_u);
        display(new_name);
        tt = tt + 1;
    end
end
%%
for tt = 1:100
    display(['C:\isbe\nailfold\models\vessel\apex_location\test_rf_' zerostr(tt, 4) '.bfs']);
end
    