This zip file includes 3 example functions to run a batch job on the cluster. These are

1)  my_function_shell_script.sh - this is the script you submit to hydra using
    qsub. For example:

        qsub -N job_name -t m-n:i path_to_script/my_function_shell_script.sh

    This will create jobs with indices [m : i : n] (using Matlab vector notation).
    So for example the command:

        qsub -N my_job -t 1-5:2 path_to_script/my_function_shell_script.sh

    Will create jobs on hydra, my_job.1, my_job.3, my_job.5 (note each of the jobs
    will also have the job ID assigned by hydra, followed by the index, for example
    843562.1, 843562.3, 843562.5)

    Within this script, the index is available for use as the variable $SGE_IDX
    The script passes the variable as an argument to my_function_wrapper (see below),
    and also uses the variable to create a separate output log for each job

2) my_function_wrapper.m - this is a Matlab function that is called by the shell script,
    creating a new instance of Matlab on Hydra.
    
    It takes as input an argument sge_index. When called by the shell script on Hydra,
    this will be the batch job index $SGE_INDEX described above.

    This function uses the index to create a different set of input conditions
    for the main function my_function (see below), for each batch job.

    In the example given, the index is used to create a different file_path in each job, so that 
    my_function is applied to a different dataset in each job.

    The function wrapper is also a useful place to set up any other input arguments
    for the main function. In this example, one of the input parameters is
    set depending on the index, thus generating a different set of
    parameters to be applied to the odd or even indexed datasets.

3) my_function.m - this is the main matlab function that actually does something
    useful. For example, loading the data, and fitting a model to it.
    