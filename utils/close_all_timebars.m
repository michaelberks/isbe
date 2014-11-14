function close_all_timebars()

% Get handles to all figures
figHandles = findall(0,'Type','figure');

for i = 1:length(figHandles)
    if strcmp(get(figHandles(i),'Tag'), 'TMWWaitbar')
        close(figHandles(i));
    end
end
