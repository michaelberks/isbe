function wavelength_data = read_wavelength_spectrum(fname)

fid = fopen(fname);
while isempty(strfind(fgetl(fid), 'Begin Processed')) && ~feof(fid)
end
if ~feof(fid); 
    wavelength_data = textscan(fid, '%f%f');
end
fclose(fid);