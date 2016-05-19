flow_results_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_results.l3\';
flow_metrics_dir = 'N:\Nailfold Capillaroscopy\Wellcome\flow_metrics\';
vessel_list = dir([flow_results_dir '*.mat']);
new_list = dir([flow_metrics_dir '*.mat']);
%%
vessel_names = {vessel_list(:).name}';
new_names = {new_list(:).name}';

missing_names = setdiff(vessel_names, new_names);
%%
missing_seqs = [];
missing_subs = [];
for i_ve = 1:length(missing_names)
    seq = missing_names{i_ve}(1:25);
    if ~ismember(seq, missing_seqs)
        missing_seqs{end+1,1} = seq; %#ok;
    end
    sub = str2double( missing_names{i_ve}(1:3) );
    if ~ismember(sub, missing_subs)
        missing_subs(end+1,1) = sub; %#ok;
    end
end