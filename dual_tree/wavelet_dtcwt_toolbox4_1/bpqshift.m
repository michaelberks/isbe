function [hnew,hnew1,hc,hcnew] = bpqshift(biort,qshift,df);
%
% Design a bandpass pair of qshift filters, with a Hilbert transform
% property and a centre frequency reduced by a factor df from a normal qshift
% high frequency filter.
% The prototype level 1 and qshift filters are specified by biort and qshift,
% and the new qshift filter hnew is the same length as the qshift filters.
% The new level 1 filter hnew1 is the same length as the highpass biort filter, h1o.
%
% hc and hcnew are returned optionally, so that the old and new impulse and freq
% responses for the complex subbands at level 2, can be plotted on fresh axes.
%
% Nick Kingsbury, Cambridge University, August 2005.

% biort = 'near_sym_b';
% qshift = 'qshift_b';
j = sqrt(-1);

load(biort)
load(qshift)

% Calculate the DTCWT highpass basis function at level 2 and plot it.
hr = [conv(h0o,upm(h1a,2)); 0]; % Add '0' to achieve 1-sample shift when hr is reversed to give hi.
hi = hr(end:-1:1);  % Imag part is real part reversed in time.
hc = hr - j*hi;  % Complex basis function.
figure(1); plot([real(hc) imag(hc) abs(hc)]); grid on

% Apply a frequency shift using a constant phase rotation rate of df radians per sample.
N = length(hc);
theta = ([1:N].' - (N+1)/2) * df;
hc2 = hc .* exp(-j*theta); % Frequency shifted basis function.

f = [-128:127]*4/256;
figure(2); plot(f,abs(fftshift(fft([hc hc2],256)))); grid on % Plot freq responses.

figure(3); plot([real(hc2) imag(hc2) abs(hc2)]); grid on % Plot desired new basis func.

% Deconvolve the level 1 lowpass filter from the real part of hc2.
hr2 = real(hc2);
hr2 = hr2(1:(end-1));  

% Form a convolution matrix of h0o.
z = zeros(length(hr2)-length(h0o),1);
Cm = toeplitz([h0o(1); z],[h0o(:); z]).';
% Strip out alternate columns of Cm to allow for upsampling of hnew.
Cm = Cm(:,1:2:end);

% Solve for hnew which gives LMS error for Cm * hnew = hr2
hnew = Cm \ hr2;
hnew = hnew - mean(hnew); % Ensure that gain = 0 at zero-freq.
% Normalise so energy of hnew is same as of h1a.
hnew = hnew * sqrt((h1a'*h1a)/(hnew'*hnew));

% Calculate the error energy.
Berr = hr2 - Cm * hnew;
err_engy = real(Berr' * Berr) 

figure(4); plot([h1a hnew]);

% Calculate new overall freq response.

hrnew = [conv(h0o,upm(hnew,2)); 0]; % Add '0' to achieve 1-sample shift when hr is reversed to give hi.
hinew = hrnew(end:-1:1);
hcnew = hrnew - j*hinew;
figure(5); plot([real(hcnew) imag(hcnew) abs(hcnew)]); grid on
figure(6); plot(f,abs(fftshift(fft([hc hcnew],256)))); grid on % Plot freq responses.

% Generate a new shifted filter for level 1.
% To be compatible with dtwavexfm, the filter must be symmetric, and create an
% approximately 1-sided freq response when combined with j times a 1-sample shift of itself.
hh = [h1o(10:19); zeros(256-19,1); h1o(1:9)]; % Zero-phase version of h1o.
HH = real(fft(hh));

% Shift HH  response down in frequency.
sh = round(df * 85); 
HHsh = [HH((sh+1):128); ones(2*sh,1)*HH(128); HH(129:(256-sh))];
HHsh = HHsh-HHsh(1); % Ensure gain = 0 at zero freq.
figure; plot(HHsh); grid on

% Apply a small amount of higher frequency attenuation
a = 0.1;
hhsh = conv(fftshift(real(ifft(HHsh))),[a 1-2*a a]);

% Select the 19-tap zero-phase filter.
hhsh1 = hhsh(130 + [-9:9]);

% Normalise so energy of hhsh1 is same as of h1o.
hhsh1 = hhsh1 * sqrt((h1o'*h1o)/(hhsh1'*hhsh1));

figure; plot([h1o hhsh1]); grid on
hcsh1=[hhsh1;0] + j*[0;hhsh1];
hc1o = [h1o;0] + j*[0;h1o];
figure; plot(abs([hc1o hcsh1])); grid on
figure; plot(f,abs(fftshift(fft([hc1o hcsh1],256)))); grid on

hnew1 = hhsh1;

return



