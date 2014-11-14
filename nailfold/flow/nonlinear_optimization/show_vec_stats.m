function show_vec_stats(vec, varname, logfile, prefix)

if ~exist('logfile','var') || isempty(logfile), logfile = 1; end
if ~exist('prefix','var') || isempty(prefix), prefix = ''; end

vec = vec(~isnan(vec));

fprintf(logfile, '%s%s:\n', prefix, varname );

prefix = [prefix '  '];
fprintf(logfile, '%sMean = %f\n', prefix, mean(vec) );
fprintf(logfile, '%sStdDev = %f\n', prefix, std(vec) );
fprintf(logfile, '%sVariance = %f\n', prefix, var(vec) );
fprintf(logfile, '%sMeanAbs = %f\n', prefix, mean(abs(vec)) );
fprintf(logfile, '%sMeanSquared = %f\n', prefix, mean(vec.^2) );
fprintf(logfile, '%sRootMeanSquared = %f\n', prefix, sqrt(mean(vec.^2)) );
