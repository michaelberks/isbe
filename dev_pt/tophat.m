function f = tophat(halfsz, halfw)

f_debug = (nargin==0 && nargout==0);
if f_debug
    halfsz = 25;
    halfw = rand*halfsz;
end

rng = -halfsz:halfsz;
f = zeros(size(rng));
for j = 1:length(rng)
    f(j) = min(1, max(0, halfw-abs(rng(j))));
end

if f_debug
    figure(1); clf;
    plot(rng,f);
    clear f;
end
