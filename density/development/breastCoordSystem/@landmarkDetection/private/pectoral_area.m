function [ area, pectoral_mask, outside_ROI ] = pectoral_area(ROI,line_rho,line_theta)
% PECTORAL_AREA returns the pectoral area (The number of pixels of the ROI
% above the line with polar coordinates line_rho and line_theta), the mask
% of this area, and whether the pectoral goes outside the ROI and thus is
% misplaced.
%
% Andreas Eilschou, 2010

pectoral_mask = logical(ROI);
outside_ROI = false;
max_zeros = 2;
number_zeros = 0;

% % Run along the pectoral line
% [N,M] = size(mask);
% for line_rho = lineprm(1):sqrt(N^2+M^2)
%     for I = N:-1:1
%         J = round((line_rho - (I-0.5) * sin(lineprm(2))) / ...
%             cos(lineprm(2)) + 0.5);
%         % original Sami: (J < 1)   -    Alternatively JM: (J < CG(1)) || (J > CG(1)+d)
%         if (J < 1) 
%             continue;
%         elseif J > M
%             break;
%         elseif mask(I,J) == 1
%             nipple = [I,J];
%             break;
%         end            
%     end
% end
%     
[N,M] = size(pectoral_mask);
for J = 1:M
    I = round((line_rho - (J-0.5) * cos(line_theta)) / ...
        sin(line_theta) + 0.5);
    if I > N
        continue;
    elseif I > 0
        % Determine whether pectoral line goes outside the ROI
        if ~outside_ROI && ROI(I,J) == 0
            number_zeros = number_zeros + 1;
            if number_zeros > max_zeros
                outside_ROI = true;
            end
        end
    else
        I = 1;
    end
    
    % Limit the mask to the pectoral area
    for l = I:N
        pectoral_mask(l,J) = false;
    end
end

% [N,M] = size(pectoral_mask);
% for J = 1:M
%     I = round((line_rho - (J-0.5) * cos(line_theta)) / ...
%         sin(line_theta) + 0.5);
%     if I > N
%         continue;
%     elseif I > 0
%         % Determine whether pectoral line goes outside the ROI
%         if ~outside_ROI && ROI(I,J) == 0
%             number_zeros = number_zeros + 1;
%             if number_zeros > max_zeros
%                 outside_ROI = true;
%             end
%         end
%     else
%         I = 1;
%     end
%     
%     % Limit the mask to the pectoral area
%     for l = I:N
%         pectoral_mask(l,J) = false;
%     end
% end
area = sum(pectoral_mask(:));

end