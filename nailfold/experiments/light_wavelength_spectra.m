data_dir = 'C:\isbe\nailfold\hardware\lightsource\NCM light source spectra\';

%Plot halogen light, green filtered halogen light and LED together
spectrum_halogen_w = read_wavelength_spectrum([data_dir 'halogen_white_6.txt']);
spectrum_halogen_g = read_wavelength_spectrum([data_dir 'halogen_green_6.txt']);
spectrum_led_g = read_wavelength_spectrum([data_dir 'LED_.txt']);
%%
figure; hold all;
plot(spectrum_halogen_w{1}, spectrum_halogen_w{2});
plot(spectrum_halogen_g{1}, spectrum_halogen_g{2});
plot(spectrum_led_g{1}, spectrum_led_g{2});
legend({'Halogen unfiltered', 'Halogen + green filter', 'Green LED'});
xlabel('Wavelength (nm)');
title('Wavelength spectra of Halogen and LED light sources');
exportfig([data_dir 'comparison.png']);
%
%Plot effect of changing temperature on Halogen light source
figure; hold all;
for ii = 1:6
    spectrum_halogen_wi = read_wavelength_spectrum([data_dir 'halogen_white_' num2str(ii) '.txt']);
    plot(spectrum_halogen_wi{1}, spectrum_halogen_wi{2});
end
xlabel('Wavelength (nm)');
title('Effect of changing temperature on unfiltered Halogen light spectrum');
exportfig([data_dir 'temperature_white.png']);

figure; hold all;
for ii = 1:6
    spectrum_halogen_wi = read_wavelength_spectrum([data_dir 'halogen_green_' num2str(ii) '.txt']);
    plot(spectrum_halogen_wi{1}, spectrum_halogen_wi{2});
end
xlabel('Wavelength (nm)');
title('Effect of changing temperature on green filtered Halogen light spectrum');
exportfig([data_dir 'temperature_green.png']);