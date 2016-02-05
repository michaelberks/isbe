function l=LineConv(pectoral)

Point1=pectoral(1,:);
Point2=pectoral(2,:);

l=cross([Point1(2),Point1(1),1],[Point2(2),Point2(1),1]);

end