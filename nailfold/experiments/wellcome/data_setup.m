%Copying over the mosaics from the network drive into a single folder on my
%local machine
local_dir = 'C:\isbe\nailfold\data\wellcome_study\full_mosaics\';
study_dir = 'N:\Nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\';
subject_dirs = dir([study_dir '0*']);

num_subs = length(subject_dirs);
recruitment_dates = zeros(num_subs,1);

for i_sub = 1:num_subs; 
    subject_dir = subject_dirs(i_sub).name;
    session_dirs = dir([study_dir subject_dir]);
    
    num_sessions = length(session_dirs);
    
    for i_sess = 1:num_sessions
        session_dir = [study_dir subject_dir '\' session_dirs(i_sess).name];
        sequences = [...
            dir([session_dir '\L*']);
            dir([session_dir '\R*'])];
        
        num_sequences = length(sequences);
        
        for i_seq = 1:num_sequences
            mosaic_name = [session_dir '\' sequences(i_seq).name '\full_mosaic.png'];
            if exist(mosaic_name, 'file');
                new_name = [local_dir ...
                    subject_dir '_' sequences(i_seq).name '_mosaic.png'];
                
                if ~exist(new_name, 'file')              
                    display(['Copying: ' mosaic_name ' to ' new_name]);
                    copyfile(mosaic_name, new_name);
                end
            end
        end
    end
end;
%%