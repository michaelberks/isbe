function [ region_map ] = region_definitions(image,region_selection,large_image)
%REGION_DEFINITIONS : Depending on the size of the image and which set of regions required; outputs the marker regions. Each region is defined by
%					  a number to denote which specific region it is. For example a small image will have 6 regions, 3 on each side. Each of the 3
%					  regions on a side will be an area of the number 1, 2 or 3 with the corresponding area on the other side of the image
%					  being a region of the same number. This is done so that markers found in a region can be identified as being found in a 
% 					  specific region so that the discrimination is made easier. Any non-region is set to zero.
% 
%
% Inputs:
%			image				image that you want to find markers in [2D array]
%			region_selection	defines which regions to output for small images [#(1 or 2)]
%
% Outputs:
%			region_map			Map of the marker regions in the image. Marker regions are areas of one's with the rest of the image being zeros
%								[2D array]
%
% Example:
%
% Notes:	Large image regions already have an added 50 pixel padding on all edges and so are written differently to the small regions in which 
%			the padding is added in this script. Hence the difference in the code for the two. 
%
%			
% See also:
%
% Created: 18-06-2013
% Author: Euan Allen
% Email : euan.allen@gmail.com

padding = 50; %how much padding to add around regions
region_map = zeros(size(image));
[Y X] = size(image);
if region_selection ==1
    if large_image %LARGE
        m_regionsL.top4 = [1037,1681,6175,6682];
        m_regionsL.top3 = [2090,2644,6175,6682];
        m_regionsL.top2 = [3105,3620,6175,6682];
        m_regionsL.top1 = [4116,4620,6175,6682];
        m_regionsL.bottom4 = [1037,1681,4,612];
        m_regionsL.bottom3 = [2090,2644,4,612];
        m_regionsL.bottom2 = [3105,3620,4,612];
        m_regionsL.bottom1 = [4116,4620,4,612];
        %Generate map from user define regions above
        region_map(m_regionsL.top1(3):m_regionsL.top1(4),m_regionsL.top1(1):m_regionsL.top1(2))=1; %Top regions
        region_map(m_regionsL.top2(3):m_regionsL.top2(4),m_regionsL.top2(1):m_regionsL.top2(2))=2;
        region_map(m_regionsL.top3(3):m_regionsL.top3(4),m_regionsL.top3(1):m_regionsL.top3(2))=3;
        region_map(m_regionsL.top4(3):m_regionsL.top4(4),m_regionsL.top4(1):m_regionsL.top4(2))=4;
        region_map(m_regionsL.bottom1(3):m_regionsL.bottom1(4),m_regionsL.bottom1(1):m_regionsL.bottom1(2))=1; %Bottom regions
        region_map(m_regionsL.bottom2(3):m_regionsL.bottom2(4),m_regionsL.bottom2(1):m_regionsL.bottom2(2))=2;
        region_map(m_regionsL.bottom3(3):m_regionsL.bottom3(4),m_regionsL.bottom3(1):m_regionsL.bottom3(2))=3;
        region_map(m_regionsL.bottom4(3):m_regionsL.bottom4(4),m_regionsL.bottom4(1):m_regionsL.bottom4(2))=4;
    else %SMALL REGION SET 1
        m_regionsS.left3 = [14-padding,373+padding,3028-padding,3236+padding];
        m_regionsS.left3(m_regionsS.left3 <1) = 1; %Does padding make region go outside of image size if so set to limit.
        if m_regionsS.left3(2)>X;  m_regionsS.left3(2) = X; end 
        if m_regionsS.left3(4)>X;  m_regionsS.left3(4) = Y; end
        
        m_regionsS.left2 = [8-padding,403+padding,2027-padding,2213+padding];
        m_regionsS.left2(m_regionsS.left2 <1) = 1;
        if m_regionsS.left2(2)>X;  m_regionsS.left2(2) = X; end
        if m_regionsS.left2(4)>X;  m_regionsS.left2(4) = Y; end
        
        m_regionsS.left1 = [32-padding,420+padding,970-padding,1264+padding];
        m_regionsS.left1(m_regionsS.left1 <1) = 1;
        if m_regionsS.left1(2)>X;  m_regionsS.left1(2) = X; end
        if m_regionsS.left1(4)>X;  m_regionsS.left1(4) = Y; end
        
        m_regionsS.right3 = [5019-padding,5217+padding,2988-padding,3227+padding];
        m_regionsS.right3(m_regionsS.right3 <1) = 1;
        if m_regionsS.right3(2)>X;  m_regionsS.right3(2) = X; end
        if m_regionsS.right3(4)>X;  m_regionsS.right3(4) = Y; end
        
        m_regionsS.right2 = [4970-padding,5181+padding,1991-padding,2264+padding];
        m_regionsS.right2(m_regionsS.right2 <1) = 1;
        if m_regionsS.right2(2)>X;  m_regionsS.right2(2) = X; end
        if m_regionsS.right2(4)>X;  m_regionsS.right2(4) = Y; end
        
        m_regionsS.right1 = [4959-padding,5152+padding,993-padding,1312+padding];
        m_regionsS.right1(m_regionsS.right1 <1) = 1;
        if m_regionsS.right1(2)>X;  m_regionsS.right1(2) = X; end
        if m_regionsS.right1(4)>X;  m_regionsS.right1(4) = Y; end
        
        region_map(m_regionsS.left1(3):m_regionsS.left1(4),m_regionsS.left1(1):m_regionsS.left1(2))=1; %left regions
        region_map(m_regionsS.left2(3):m_regionsS.left2(4),m_regionsS.left2(1):m_regionsS.left2(2))=2;
        region_map(m_regionsS.left3(3):m_regionsS.left3(4),m_regionsS.left3(1):m_regionsS.left3(2))=3;
        region_map(m_regionsS.right1(3):m_regionsS.right1(4),m_regionsS.right1(1):m_regionsS.right1(2))=1; %right regions
        region_map(m_regionsS.right2(3):m_regionsS.right2(4),m_regionsS.right2(1):m_regionsS.right2(2))=2;
        region_map(m_regionsS.right3(3):m_regionsS.right3(4),m_regionsS.right3(1):m_regionsS.right3(2))=3;
    end
elseif region_selection ==2
    if length(image)>6000 %LARGE
        m_regionsL.top4 = [1037,1681,6175,6682];
        m_regionsL.top3 = [2090,2644,6175,6682];
        m_regionsL.top2 = [3105,3620,6175,6682];
        m_regionsL.top1 = [4116,4620,6175,6682];
        m_regionsL.bottom4 = [1037,1681,4,612];
        m_regionsL.bottom3 = [2090,2644,4,612];
        m_regionsL.bottom2 = [3105,3620,4,612];
        m_regionsL.bottom1 = [4116,4620,4,612];
        %Generate map from user define regions above
        region_map(m_regionsL.top1(3):m_regionsL.top1(4),m_regionsL.top1(1):m_regionsL.top1(2))=1;
        region_map(m_regionsL.top2(3):m_regionsL.top2(4),m_regionsL.top2(1):m_regionsL.top2(2))=2;
        region_map(m_regionsL.top3(3):m_regionsL.top3(4),m_regionsL.top3(1):m_regionsL.top3(2))=3;
        region_map(m_regionsL.top4(3):m_regionsL.top4(4),m_regionsL.top4(1):m_regionsL.top4(2))=4;
        region_map(m_regionsL.bottom1(3):m_regionsL.bottom1(4),m_regionsL.bottom1(1):m_regionsL.bottom1(2))=1;
        region_map(m_regionsL.bottom2(3):m_regionsL.bottom2(4),m_regionsL.bottom2(1):m_regionsL.bottom2(2))=2;
        region_map(m_regionsL.bottom3(3):m_regionsL.bottom3(4),m_regionsL.bottom3(1):m_regionsL.bottom3(2))=3;
        region_map(m_regionsL.bottom4(3):m_regionsL.bottom4(4),m_regionsL.bottom4(1):m_regionsL.bottom4(2))=4;
    else %SMALL REGION SET 2
        m_regionsS.left3 = [91-padding,416+padding,2404-padding,2603+padding];
        m_regionsS.left3(m_regionsS.left3 <1) = 1;
        if m_regionsS.left3(2)>X;  m_regionsS.left3(2) = X; end
        if m_regionsS.left3(4)>X;  m_regionsS.left3(4) = Y; end
        
        m_regionsS.left2 = [115-padding,418+padding,1395-padding,1626+padding];
        m_regionsS.left2(m_regionsS.left2 <1) = 1;
        if m_regionsS.left2(2)>X;  m_regionsS.left2(2) = X; end
        if m_regionsS.left2(4)>X;  m_regionsS.left2(4) = Y; end
        
        m_regionsS.left1 = [199-padding,419+padding,423-padding,682+padding];
        m_regionsS.left1(m_regionsS.left1 <1) = 1;
        if m_regionsS.left1(2)>X;  m_regionsS.left1(2) = X; end
        if m_regionsS.left1(4)>X;  m_regionsS.left1(4) = Y; end
        
        m_regionsS.right3 = [4982-padding,5303+padding,2360-padding,2609+padding];
        m_regionsS.right3(m_regionsS.right3 <1) = 1;
        if m_regionsS.right3(2)>X;  m_regionsS.right3(2) = X; end
        if m_regionsS.right3(4)>X;  m_regionsS.right3(4) = Y; end
        
        m_regionsS.right2 = [4974-padding,5257+padding,1313-padding,1619+padding];
        m_regionsS.right2(m_regionsS.right2 <1) = 1;
        if m_regionsS.right2(2)>X;  m_regionsS.right2(2) = X; end
        if m_regionsS.right2(4)>X;  m_regionsS.right2(4) = Y; end
        
        m_regionsS.right1 = [4967-padding,5204+padding,459-padding,689+padding];
        m_regionsS.right1(m_regionsS.right1 <1) = 1;
        if m_regionsS.right1(2)>X;  m_regionsS.right1(2) = X; end
        if m_regionsS.right1(4)>X;  m_regionsS.right1(4) = Y; end
        
        region_map(m_regionsS.left1(3):m_regionsS.left1(4),m_regionsS.left1(1):m_regionsS.left1(2))=1;
        region_map(m_regionsS.left2(3):m_regionsS.left2(4),m_regionsS.left2(1):m_regionsS.left2(2))=2;
        region_map(m_regionsS.left3(3):m_regionsS.left3(4),m_regionsS.left3(1):m_regionsS.left3(2))=3;
        region_map(m_regionsS.right1(3):m_regionsS.right1(4),m_regionsS.right1(1):m_regionsS.right1(2))=1;
        region_map(m_regionsS.right2(3):m_regionsS.right2(4),m_regionsS.right2(1):m_regionsS.right2(2))=2;
        region_map(m_regionsS.right3(3):m_regionsS.right3(4),m_regionsS.right3(1):m_regionsS.right3(2))=3;
    end
end
end