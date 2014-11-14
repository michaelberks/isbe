function smoothed = smooth(input,width)

dim = size(input);

smoothed=zeros(1,dim(2));

% pad out array with zeros
temp = zeros(1,dim(2)+2*width);
temp(width+1:dim(2)+width)=input;

for i=width+1:dim(2)+width
    smoothed(i-width)=mean(temp(i-round(width/2):i+round(width/2)));
end
