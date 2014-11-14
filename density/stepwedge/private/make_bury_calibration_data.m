function calibration_data = make_bury_calibration_data(save_path)

calibration_data(1).M_mag = [8.35 8.35 8.36]; %Note we now miss off the first marker in 18x24 films: 8.10 
calibration_data(2).M_mag = [10.32 10.33 10.41 10.42];
calibration_data(1).C_mag = [4472.20 4477.70 4479.50]; %Note we now miss off the first marker in 18x24 films: 4475.20
calibration_data(2).C_mag = [5515.30 5517.50 5522.60 5527.50];
calibration_data(1).stepwedge = zeros(1,3);

gt{1} = [0 5 10 15 20];
gt{2} = [0 10 15 20 25 30];
gt{3} = [10 20 25 30 35];
gt{4} = [0 10 15 20 25 30 35 40 50];
gt{5} = [0 10 15 20 35 40 45 50];
gt{6} = [0 5 10 15 20 25 30 35 40 45];
gt{7} = [0 10 20 30];
gt{8} = [0 10 20 30];
gt{9} = [0 10 20];
gt{10} = [0 10 20];

% y-coordinates = stepwedge thickness

sh{1} = [0.81 0.97 1.14 1.34 1.51];
sh{2} = [1.29 1.67 1.87 2.06 2.30 2.50];
sh{3} = [2.24 2.75 3.01 3.31 3.67];
sh{4} = [2.43 2.95 3.33 3.71 3.90 4.21 4.82 5.21 6.86];
sh{5} = [3.29 3.78 4.18 4.60 6.16 6.70 7.59 9.07];
sh{6} = [3.97 4.38 4.83 5.35 6.02 6.66 7.54 8.63 10.60 13.99];
sh{7} = [4.93 5.62 7.26 9.77];
sh{8} = [5.90 6.79 8.84 12.25];
sh{9} = [6.95 8.07 10.61];
sh{10} = [8.07 9.48 12.56];

row = 1;
for ii = 1:10
    bt = 10*(ii + 1);
    for jj = 1:length(gt{ii})   
        calibration_data(1).stepwedge(row,:) = [bt gt{ii}(jj) sh{ii}(jj)];
        row = row+1;
    end
end

calibration_data(2).stepwedge = calibration_data(1).stepwedge;

if nargin > 0
    save(save_path, 'calibration_data');
end