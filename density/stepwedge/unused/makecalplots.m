% results of calibration investigations
% results are area % area>5mm % vol % vol cm3, vol_b cm3
% volumes need multiplying by 25 to get to cm3
old=[68.92 39.97   10.5 6.11     58.19
66.92 40.51   8.89 6.26     70.41
64.95 34.94   5.95 3.36     56.47
72.82 55.08  12.22 5.31     43.45
69.26 21.83   4.18 4.32    103.3
79.48 55.05  18.46 2.55    13.81
75.61 69.04  31.45 8.54    27.15  
72.16 63.57  26.69 3.47    13.00  ];
new=[67.51 35.81   8.84 5.11    58.19
62.94 37.54   8.31 6.26    75.33
62.82 28.52   5.07 2.84    56.47
73.04 55.02  12.30 5.38    43.45
65.48 16.71   3.42 3.52    103.3
78.09 58.95  19.62 2.77    13.81
72.74 65.00  27.35 7.54    27.15
70.71 62.34  29.34 3.40    11.59];
newdig=[82.62 74.37  19.67 10.88    55.32 
73.37 49.75  11.18 7.80     69.78
70.45 51.58   8.44 4.47     53.02    
79.40 68.27  17.72 7.31     41.25
73.93 23.80   4.60 4.59     99.78
86.28 79.05  30.83 3.85     12.52
80.12 72.52  33.81 8.77     25.93   
86.00 78.76  50.13 4.70     9.37];
newcal=[80.28  70.83 18.41 10.57    57.41
72.49 42.32  10.04 7.15     71.24
69.72 51.74   8.54 4.66     54.59
71.86 51.11  11.20 4.59     40.98
68.89 23.57   4.74 4.84     102.11
76.38 43.03  15.15 2.07     13.65
83.83 74.57  34.47 8.70     25.24
85.31 77.53  49.81 4.76     9.57];
cal2004nosw=[77.30  64.60 15.95 9.40     58.93
65.17 33.50   8.27 5.93     71.70 
62.71 28.00   4.73 2.70     57.11
75.82 61.29  12.59 5.44     43.24
52.95 11.42   2.47 2.57    104.12
76.38 43.03  15.15 2.07    13.65
77.40 67.37  26.63 4.17    26.39
63.48 49.56  20.63 2.63   12.74];
cal2004=[79.09  67.08 16.30 9.51     55.8
68.47 35.69   8.35 5.98     71.74
65.12 32.55   5.29 2.99     56.60
76.56 60.94  12.53 5.44     43.42
58.21 12.98   2.84 2.94    103.69
79.10 47.35  15.88 2.19    13.76
78.53 68.49  26.26 6.90    26.28
65.13 50.83  20.27 2.57   12.66 ];
tentimes=[85.81             28.97          12.71          3.68            77.55
87.42             32.47          12.45          4.04            80.72
84.89             26.32          12.87          3.39            73.18
86.17             31.51          12.42          3.91            79.42
85.13             28.08          12.64          3.55            76.64
85.16             28.72          12.56          3.61            77.20
85.71             29.50          12.62          3.72            78.30
84.87             27.14          12.70          3.45            75.59
86.11             31.01          12.54          3.89            79.25
84.29             26.56          12.83          3.41            74.50
];
old(:,4)=old(:,4)*25;
old(:,5)=old(:,5)*25;
new(:,4)=new(:,4)*25;
new(:,5)=new(:,5)*25;
newdig(:,4)=newdig(:,4)*25;
newdig(:,5)=newdig(:,5)*25;
newcal(:,4)=newcal(:,4)*25;
newcal(:,5)=newcal(:,5)*25;
cal2004(:,4)=cal2004(:,4)*25;
cal2004(:,5)=cal2004(:,5)*25;
tentimes(:,4)=tentimes(:,4)*25;
tentimes(:,5)=tentimes(:,5)*25;
% get difference between my results and old results
diffnew=new-old;
diffnewdig=newdig-old;
diffnewcal=newcal-old;
diffcal2004=cal2004-old;
% calculate normalised difference so can consider all images together
diffnew=diffnew./old;
diffnewdig=diffnewdig./old;
diffnewcal=diffnewcal./old;
diffcal2004=diffcal2004./old;
% plot some bar charts - first, % glandular vol
barstatspv(:,1)=diffnew(:,3);
barstatspv(:,2)=diffnewdig(:,3);
barstatspv(:,3)=diffnewcal(:,3);
barstatspv(:,4)=diffcal2004(:,3);
figure;h=bar(barstatspv,'grouped');
%set(h,'Name','percent v_g');
% now glandular vol
barstatsgv(:,1)=diffnew(:,4);
barstatsgv(:,2)=diffnewdig(:,4);
barstatsgv(:,3)=diffnewcal(:,4);
barstatsgv(:,4)=diffcal2004(:,4);
figure;h=bar(barstatsgv,'grouped');
%set(h,'Name',' v_g cm3');
% now breast vol
barstatsbv(:,1)=diffnew(:,5);
barstatsbv(:,2)=diffnewdig(:,5);
barstatsbv(:,3)=diffnewcal(:,5);
barstatsbv(:,4)=diffcal2004(:,5);
figure;h=bar(barstatsbv,'grouped');
%set(h,'Name','v_b cm3');
figure;bar(tentimes(:,3));