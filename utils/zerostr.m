function newstr = zerostr(num, len)
    
    newstr = sprintf(['%0' num2str(len) 'u'], num);
        
end