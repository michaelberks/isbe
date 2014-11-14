[num txt raw] = xlsread('C:\isbe\density\misc\mamm Dates1.xls');
new_data = cell(length(txt),1);

for ii = 1:length(txt);
    if ~isempty(txt{ii})
        temp = txt{ii};
        odd_chars = ~isstrprop(temp, 'punct') & ~isstrprop(temp, 'digit');
        temp(odd_chars) = [];
        
        if rem(length(temp), 10) || isempty(temp)
            %display(['Cell ', num2str(ii), 'not divisible by 10: ', temp]);
        else
            ll = 1;
            for kk = 1:10:length(temp)-9
                new_data{ii, ll} = temp(kk:kk+9);
                ll = ll+1;
            end

        end
    end
end

xlswrite('C:\isbe\density\misc\mamm Dates1.xls', new_data, 2);