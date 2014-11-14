% setbpfilters.m
%
% Set up the bandpass versions of the standard filters to implement
% the bandpass 45 degree subbands.
%
% Nick Kingsbury, Cambridge University, August 2005.

biort = 'near_sym_b';
qshift = 'qshift_b';

load(biort)
load(qshift)

close all;

% Apply a downward frequency shift to the highpass subband at level 2.
df = 0.293; % Frequency shift in radians per sample at level 0.
[h2a,h2o,hc,hcnew] = bpqshift(biort,qshift,df);

h2b = h2a(end:-1:1);
g2a = h2b; g2b = h2a;
g2o = h2o;

save([qshift '_bp'],'h0a','h0b','h1a','h1b','h2a','h2b','g0a','g0b','g1a','g1b','g2a','g2b');

save([biort '_bp'],'h0o','h1o','h2o','g0o','g1o','g2o');

%showbpfilters





