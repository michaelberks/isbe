%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%        Experiment to build RF classifier using the ISBE code
%--------------------------------------------------------------------------

%1. Create some images (e.g. containing synthetic linear structures on real
% mammographic backgrounds) and copy to isbe-san1

%2. On hydra, submit a batch of jobs, each of which creates a given number
% of trees - the matlab function to do this is:
%   build_rf(job_idx, image_dir, n_samples, n_trees, d)
% This function can be called on hydra by executing the script build_rf.sh
%
% So, for example, to create 10 jobs each of which builds 50 trees, where
% each tree randomly selects 200,000 samples from the images in image_dir
% and uses a subset of 10 randomly selected dimensions and each tree split,
% edit the 3rd line of build_rf.sh to read:
%   matlab -nodisplay -nojvm -r "build_rf($SGE_TASK_ID, 'image_dir', 2e5, 10, 10)" > logs/build_rf$SGE_TASK_ID.log
%
% Now execute:
%   qsub -N rf_job_name -t 1-10 matlab_code/trunk/hydra/build_rf.sh

%3. When the jobs have finished, run mb_combine_rfs to merge the separate
% forests. You can also use this function to copy the merged (forest and
% its trees) from isbe-san1 to your local machine

%4. The forest can then be used to classify pixels in a real mammogram,
%effectively performing line detection