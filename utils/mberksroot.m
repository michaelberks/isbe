function myroot = mberksroot
% get correct root variable for pc/cluster    
    if ispc %work/home desktop
        myroot = 'C:/isbe/dev/';
    else % hydra cluster
        myroot = '/san/staff/student/mberks/dev/';
    end