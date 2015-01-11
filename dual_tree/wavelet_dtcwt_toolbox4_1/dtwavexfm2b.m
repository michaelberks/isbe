function [Yh,Yscale] = dtwavexfm2b(X,nlevels,biort,qshift,quiet)

% Function to perform a n-level DTCWT-2D decompostion on a 2D matrix X
%
% [Yl,Yh,Yscale] = dtwavexfm2b(X,nlevels,biort,qshift);
%
%     X -> 2D real matrix/Image
%
%     nlevels -> No. of levels of wavelet decomposition
%
%     biort ->  'antonini'   => Antonini 9,7 tap filters.
%               'legall'     => LeGall 5,3 tap filters.
%               'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               'near_sym_b' => Near-Symmetric 13,19 tap filters.
%               'near_sym_b_bp' => Near-Symmetric 13,19 tap filters with extra bandpass filters.
%
%     qshift -> 'qshift_06' => Quarter Sample Shift Orthogonal (Q-Shift) 10,10 tap filters, 
%                              (only 6,6 non-zero taps).
%               'qshift_a' =>  Q-shift 10,10 tap filters,
%                              (with 10,10 non-zero taps, unlike qshift_06).
%               'qshift_b' => Q-Shift 14,14 tap filters.
%               'qshift_c' => Q-Shift 16,16 tap filters.
%               'qshift_d' => Q-Shift 18,18 tap filters.
%               'qshift_b_bp' => Q-Shift 14,14 tap filters with extra bandpass filters.
%               
%
%     Yl     -> The real lowpass image from the final level
%     Yh     -> A cell array containing the 6 complex highpass subimages for each level.
%     Yscale -> This is an OPTIONAL output argument, that is a cell array containing 
%               real lowpass coefficients for every scale.
%
% 
% Example: [Yl,Yh] = dtwavexfm2(X,3,'near_sym_b','qshift_b');
% performs a 3-level transform on the real image X using the 13,19-tap filters 
% for level 1 and the Q-shift 14-tap filters for levels >= 2.
% Version 2b allows for bandpass filters for 45 degree subbands to give
% improved rotational symmetry during analysis (but destroys P-R property).
%
% Nick Kingsbury, Cambridge University, August 2005

if ~isa(X, 'double')
    X = double(X);
end

%set default filters
if nargin < 3
    biort = 'near_sym_b_bp';
    qshift = 'qshift_b_bp';
end
if nargin < 5
    quiet = true;
end

if isstr(biort) && isstr(qshift)		%Check if the inputs are strings
   biort_exist = exist([biort '.mat'], 'file');
   qshift_exist = exist([qshift '.mat'], 'file');
   if 0%biort_exist && qshift_exist;        		%Check to see if the inputs exist as .mat files
      load(biort);
      load(qshift);
   else
      [h0o h1o g0o g1o h2o g2o h0a g0a h1a h0b h1b g1a g0b g1b h2a h2b g2a g2b] = make_default_filts();
      %error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEXFM2 for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEXFM2.');
end 

orginal_size = size(X);

if ndims(X) >= 3;
   error(sprintf('The entered image is %dx%dx%d, please enter each image slice separately.',orginal_size(1),orginal_size(2),orginal_size(3)));
end

bp_qsh = (exist('h2a')==1) && (exist('h2b')==1);
bp_lev1 = (exist('h2o')==1); 
hi_lev1 = (exist('h1o')==1);

% The next few lines of code check to see if the image is odd in size, if so an extra ...
% row/column will be added to the bottom/right of the image
initial_row_extend = 0;  %initialise
initial_col_extend = 0;
if any(rem(orginal_size(1),2)), %if sx(1) is not divisable by 2 then we need to extend X by adding a row at the bottom
   X = [X; X(end,:)];           %Any further extension will be done in due course.
   initial_row_extend = 1;
end
if any(rem(orginal_size(2),2)), 	%if sx(2) is not divisable by 2 then we need to extend X by adding a col to the left
   X = [X X(:,end)];          %Any further extension will be done in due course.
   initial_col_extend = 1;
end
extended_size = size(X);

if nlevels == 0, return; end

%initialise
Yh=cell(nlevels+1,1);
if nargout == 3
   Yscale=cell(nlevels,1);   %this is only required if the user specifies a third output component.
end

S = [];
sx = size(X);
if nlevels >= 1,
   
   % Do odd top-level filters on cols.
   Lo = colfilter(X,h0o).';
   if hi_lev1, Hi = colfilter(X,h1o).'; end
   if bp_lev1, Ba = colfilter(X,h2o).'; end
   
   % Do odd top-level filters on rows.
   LoLo = colfilter(Lo,h0o).';			% LoLo
   if hi_lev1,
       Yh{1} = zeros([size(LoLo)/2  6]);
       Yh{1}(:,:,[1 6]) = q2c(colfilter(Hi,h0o).');			% Horizontal pair
       Yh{1}(:,:,[3 4]) = q2c(colfilter(Lo,h1o).');			% Vertical pair
       if bp_lev1,
           Yh{1}(:,:,[2 5]) = q2c(colfilter(Ba,h2o).');	    % Diagonal bandpass pair
       else
           Yh{1}(:,:,[2 5]) = q2c(colfilter(Hi,h1o).');	    % Diagonal highpass pair
       end
   end
   S = [ size(LoLo) ;S];
   if nargout == 3
      Yscale{1} = LoLo;
   end
end

if nlevels >= 2;
   for level = 2:nlevels;
      [row_size col_size] = size(LoLo);
      if any(rem(row_size,4)),		% Extend by 2 rows if no. of rows of LoLo are not divisable by 4;
         LoLo = [LoLo(1,:); LoLo; LoLo(end,:)];
      end 
      if any(rem(col_size,4)),		% Extend by 2 cols if no. of cols of LoLo are not divisable by 4;
         LoLo = [LoLo(:,1)  LoLo  LoLo(:,end)];
      end 
      
      % Do even Qshift filters on rows.
      Lo = coldfilt(LoLo,h0b,h0a).';
      Hi = coldfilt(LoLo,h1b,h1a).';
      if bp_qsh, Ba = coldfilt(LoLo,h2b,h2a).'; end
      
      % Do even Qshift filters on columns.
      LoLo = coldfilt(Lo,h0b,h0a).';	%LoLo
      LoLo = LoLo / 2; %account for factoring in decimation
      
      Yh{level} = zeros([size(LoLo)/2  6]);
      Yh{level}(:,:,[1 6]) = q2c(coldfilt(Hi,h0b,h0a).');	% Horizontal
      Yh{level}(:,:,[3 4]) = q2c(coldfilt(Lo,h1b,h1a).');	% Vertical
      if bp_qsh,
          Yh{level}(:,:,[2 5]) = q2c(coldfilt(Ba,h2b,h2a).');	% Diagonal bandpass 
      else
          Yh{level}(:,:,[2 5]) = q2c(coldfilt(Hi,h1b,h1a).');	% Diagonal highpass
      end
      S = [ size(LoLo) ;S];
      if nargout == 3
         Yscale{level} = LoLo;
      end
   end
end

%Save final low-pass band in output
Yh{end} = LoLo;

if ~quiet
    if initial_row_extend == 1 && initial_col_extend == 1;
       warning(sprintf(' \r\r The image entered is now a %dx%d NOT a %dx%d \r The bottom row and rightmost column have been duplicated, prior to decomposition. \r\r ',...
          extended_size(1),extended_size(2),orginal_size(1),orginal_size(2)));
    elseif initial_row_extend == 1 ;
       warning(sprintf(' \r\r The image entered is now a %dx%d NOT a %dx%d \r Row number %d has been duplicated, and added to the bottom of the image, prior to decomposition. \r\r',...
          extended_size(1),extended_size(2),orginal_size(1),orginal_size(2),orginal_size(1)));
    elseif initial_col_extend == 1;
       warning(sprintf(' \r\r The image entered is now a %dx%d NOT a %dx%d \r Col number %d has been duplicated, and added to the right of the image, prior to decomposition. \r\r',...
          extended_size(1),extended_size(2),orginal_size(1),orginal_size(2),orginal_size(2)));
    end
end
return

%==========================================================================================
%						**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function z = q2c(y)

% function z = q2c(y)
% Convert from quads in y to complex numbers in z.

sy = size(y);
t1 = 1:2:sy(1); t2 = 1:2:sy(2);
j2 = sqrt([0.5 -0.5]);

% Arrange pixels from the corners of the quads into
% 2 subimages of alternate real and imag pixels.
%  a----b
%  |    |
%  |    |
%  c----d

% Combine (a,b) and (d,c) to form two complex subimages. 
p = y(t1,t2)*j2(1) + y(t1,t2+1)*j2(2);     % p = (a + jb) / sqrt(2)
q = y(t1+1,t2+1)*j2(1) - y(t1+1,t2)*j2(2); % q = (d - jc) / sqrt(2)

% Form the 2 subbands in z.
z = cat(3,p-q,p+q);

return;

%==========================================================================================
%						**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function [h0o h1o g0o g1o h2o g2o h0a g0a h1a h0b h1b g1a g0b g1b h2a h2b g2a g2b] = make_default_filts()

h0o = [-0.0017578; 0; 0.022266; -0.046875; -0.048242; 0.29688; 0.55547; 0.29688; -0.048242; -0.046875; 0.022266; 0; -0.0017578; ];
h1o = [-7.0626e-05; 0; 0.0013419; -0.0018834; -0.0071568; 0.023856; 0.055643; -0.051688; -0.29976; 0.55943; -0.29976; -0.051688; 0.055643; 0.023856; -0.0071568; -0.0018834; 0.0013419; 0; -7.0626e-05; ];
g0o = [7.0626e-05; 0; -0.0013419; -0.0018834; 0.0071568; 0.023856; -0.055643; -0.051688; 0.29976; 0.55943; 0.29976; -0.051688; -0.055643; 0.023856; 0.0071568; -0.0018834; -0.0013419; 0; 7.0626e-05; ];
g1o = [-0.0017578; 0; 0.022266; 0.046875; -0.048242; -0.29688; 0.55547; -0.29688; -0.048242; 0.046875; 0.022266; 0; -0.0017578; ];
h2o = [-0.00036825; -0.00062225; -7.8178e-05; 0.0041858; 0.0081918; -0.0074233; -0.061538; -0.14816; -0.11708; 0.65291; -0.11708; -0.14816; -0.061538; -0.0074233; 0.0081918; 0.0041858; -7.8178e-05; -0.00062225; -0.00036825; ];
g2o = [-0.00036825; -0.00062225; -7.8178e-05; 0.0041858; 0.0081918; -0.0074233; -0.061538; -0.14816; -0.11708; 0.65291; -0.11708; -0.14816; -0.061538; -0.0074233; 0.0081918; 0.0041858; -7.8178e-05; -0.00062225; -0.00036825; ];

h0a = [0.0032531; -0.0038832; 0.03466; -0.038873; -0.1172; 0.2753; 0.75615; 0.56881; 0.011866; -0.10671; 0.023825; 0.017025; -0.0054395; -0.0045569;];
g0a = [-0.0045569; -0.0054395; 0.017025; 0.023825; -0.10671; 0.011866; 0.56881; 0.75615; 0.2753; -0.1172; -0.038873; 0.03466; -0.0038832; 0.0032531;];
h1a = [-0.0045569; 0.0054395; 0.017025; -0.023825; -0.10671; -0.011866; 0.56881; -0.75615; 0.2753; 0.1172; -0.038873; -0.03466; -0.0038832; -0.0032531;];
h0b = [-0.0045569; -0.0054395; 0.017025; 0.023825; -0.10671; 0.011866; 0.56881; 0.75615; 0.2753; -0.1172; -0.038873; 0.03466; -0.0038832; 0.0032531;];
h1b = [-0.0032531; -0.0038832; -0.03466; -0.038873; 0.1172; 0.2753; -0.75615; 0.56881; -0.011866; -0.10671; -0.023825; 0.017025; 0.0054395; -0.0045569;];
g1a = [-0.0032531; -0.0038832; -0.03466; -0.038873; 0.1172; 0.2753; -0.75615; 0.56881; -0.011866; -0.10671; -0.023825; 0.017025; 0.0054395; -0.0045569;];
g0b = [0.0032531; -0.0038832; 0.03466; -0.038873; -0.1172; 0.2753; 0.75615; 0.56881; 0.011866; -0.10671; 0.023825; 0.017025; -0.0054395; -0.0045569;];
g1b = [-0.0045569; 0.0054395; 0.017025; -0.023825; -0.10671; -0.011866; 0.56881; -0.75615; 0.2753; 0.1172; -0.038873; -0.03466; -0.0038832; -0.0032531;];
h2a = [-2.4356e-05; -0.0095951; -0.025455; -0.026369; -0.0076247; 0.26269; 0.43679; -0.83814; -0.044765; 0.17324; 0.061445; 0.02101; -0.00043292; -0.0027717;];
h2b = [-0.0027717; -0.00043292; 0.02101; 0.061445; 0.17324; -0.044765; -0.83814; 0.43679; 0.26269; -0.0076247; -0.026369; -0.025455; -0.0095951; -2.4356e-05;];
g2a = [-0.0027717; -0.00043292; 0.02101; 0.061445; 0.17324; -0.044765; -0.83814; 0.43679; 0.26269; -0.0076247; -0.026369; -0.025455; -0.0095951; -2.4356e-05;];
g2b = [-2.4356e-05; -0.0095951; -0.025455; -0.026369; -0.0076247; 0.26269; 0.43679; -0.83814; -0.044765; 0.17324; 0.061445; 0.02101; -0.00043292; -0.0027717;];
