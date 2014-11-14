function [output]=CalculateInteractionFunction(input)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% function CalculateInteractionFunction
%
% syntax: [output]=CalculateInteractionFunction(input)
% 
% Output: 
%       output: the output image.
%            
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

 
%create and get the size of the mask
mask=ones(5,5);      
[sz sz]= size(mask);
  
%stretch the data to avoid the annoying black window of the mask
%related operations
for i=1:(sz-1)/2      
    input=stretch_data_minusonepad(input); 
end
       
%get the new dimensions of the input
[dimy dimx]= size(input);
output=zeros(size(input));
       
% create a vector to store the calculations    
c=zeros(length(mask),1);

%main loop
for i=(1+ (sz-1)/2):dimx-(sz-1)/2                % x dimension
    for j=(1+ (sz-1)/2):dimy-(sz-1)/2            % y dimension
    counter=0;
    
        for mcr2=1:sz              % x mask dimension                         
            for mcr1=1:sz           % y mask dimension
  
            c1=(mcr1-((sz-1)/2)-1);
            c2=(mcr2-((sz-1)/2)-1);
    
            counter=counter+1;
            %intensity
%             x=abs(input(j+c2, i+c1)-input(j, i));
                   
            %euclidean distance based on true coordinates
            %sqrt((j+c2-j)^2 + (i+c1-i)^2)
            x=sqrt(c2^2+c1^2);                  
                 if x<2
                       c(counter)= (1/2)*(2-(abs(x)))+((1/pi)*sin(pi*abs(x)));                                              
                 else
                       c(counter)=0; 
                 end
            end
        end
 
        c_sum=sum(c(:));
        output(j,i)=c_sum;
    end%y dimension
end%x dimension        
    

for i=1:(sz-1)/2
    output=crop_data(output);     
end 