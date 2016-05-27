function myroot = toyotaroot(m_drive)
% get correct root variable for pc/cluster    
    if ispc %work/home desktop
        if nargin < 1 
            m_drive = 0;
        end
        if m_drive
            myroot = '';
        else
            switch get_username
                case {'mberks', 'Michael Berks', 'momeemb2'}
                    myroot = 'C:\isbe\toyota\';
                otherwise
                    myroot = 'C:\isbe\toyota\';
            end
        end
    else %
        switch get_username
            case {'mberks'} %Yasmina's PC
                myroot = '/home/mberks/isbe/toyota/';
            case 'momeemb2' %CSF
                myroot = 'scratch/toyota/';
            otherwise
                myroot = 'scratch/toyota/';
        end
        
        %myroot = '/home/mberks/asymmetry_project/';
    end