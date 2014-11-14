rad_list = dir('C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*.mat');

for ii = 1:length(rad_list)
    c_name = ['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad_list(ii).name];
    m_name = ['M:\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad_list(ii).name];
    if ~exist(m_name, 'file')
        copyfile(c_name, m_name);
    end
end
%%
rad_list = dir('Z:\asymmetry_project\data\radial_maps\2004_screening\normals\*.mat');

for ii = 1:length(rad_list)
    c_name = ['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\normals\' rad_list(ii).name];
    z_name = ['Z:\asymmetry_project\data\radial_maps\2004_screening\normals\' rad_list(ii).name];
    if ~exist(c_name, 'file')
        copyfile(z_name, m_name);
    end
end