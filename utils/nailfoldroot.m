function myroot = nailfoldroot(m_drive)
% get correct root variable for pc/cluster    
    if ispc %work/home desktop
        if nargin < 1 
            m_drive = 0;
        end
        if m_drive
            myroot = 'M:\nailfold\';
        else
            switch get_username
                case {'mberks', 'Michael Berks', 'momeemb2'}
                    myroot = 'C:\isbe\nailfold\';
                case 'ptresadern'
                    myroot = 'U:\projects\mammography\';
                otherwise
                    myroot = 'E:\nailfold\';
            end
        end
    else % hydra cluster
        switch get_username
            case {'mberks', 'ptresadern'}
                myroot = '/san/images/asym/nailfold/';
            case 'momeemb2'
                myroot = 'scratch/nailfold/';
            otherwise
                myroot = '/san/images/asym/nailfold/';
        end
        
        %myroot = '/home/mberks/asymmetry_project/';
    end