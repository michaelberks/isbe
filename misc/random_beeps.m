function random_beeps(num_beeps)

cmd{1} = 'powershell -inputformat none [System.Media.SystemSounds]::Beep.Play()'; 
cmd{2} = 'powershell -inputformat none [System.Media.SystemSounds]::Asterisk.Play()'; 
cmd{3} = 'powershell -inputformat none [System.Media.SystemSounds]::Hand.Play()'; 
cmd{4} = 'powershell -inputformat none [System.Media.SystemSounds]::Exclamation.Play()'; 

if ischar(num_beeps)
    num_beeps = str2double(num_beeps);
end
for i_b = 1:num_beeps
    r = ceil(4*rand);    
    [status,result] = system(cmd{r});
    pause(1);
    %display(status);
    %display(result);
end
        