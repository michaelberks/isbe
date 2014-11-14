% QSUB
% Intsructions for submitting a batch clustering job with qsub on hydra
%
% 1) Build cluster idx files using MB_SAVE_CLUSTER_IDX either on hydra or
%    offline
% 2) Make sure images, temp dirs and cluster idxs are transferred to hydra
% 3) Amend file paths in CLUSTER_WRAPPER and FINAL_CLUSTER_WRAPPER for
%    current dataset
% 4) Submit the initial divide stage jobs using the following command:
%    "qsub -N divide_name -t 1-n:1 matlab_code/trunk/cluster/divide_cluster_script.sh" 
%    where n is the number of divided groups to cluster
% 5) Submit the final clustering job, held until divide stage completes
%    "qsub -N final_name -hold_jid divide_name* matlab_code/trunk/cluster/final_cluster_script.sh"
%
% See also MB_CLUSTER_ONCE and MB_CLUSTER_FINAL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

