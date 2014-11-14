function form_name = make_TPform_pdf(full_name, second_name, nhs_id, dob, sx_id, app_date_curr, app_date_prev, form_folder)

%Notes - as of 4th November 2009, marker 4 is used in place of marker 3 to
%distinguish between copies of the form in which the order of lines has
%changed

%load in markers (marker1, marker2, marker3)
load components.mat marker*

% Build lines
line_scale = true(200, 600);
for ii = 1:4
    offset = 50*(ii-1);
    line_scale((24:26) + offset, 101:500) = 0;
end

yes_no = true(75, 600);
yes_no([25 50], 300:325) = 0;
yes_no(25:50, [300 325]) = 0;
yes_no([25 50], 400:425) = 0;
yes_no(25:50, [400 425]) = 0;

%Get binary no for nhs_id number
nhs_copy = nhs_id;
nhs_bin = zeros(1,34);
for ii = 33:-1:0
    nhs_bin(34 - ii) = 2^ii <= nhs_copy;
    nhs_copy = rem(nhs_copy, 2^ii);
end

%Convert date string into numbers dd, mm and yy
dob_dd = str2double(dob(1:2));
dob_mm = str2double(dob(4:5));
dob_yy = str2double(dob(7:end)); %assume this could also be yyyy

%Get binary no for dd
dd_bin = zeros(1,5);
for ii = 4:-1:0
    dd_bin(5 - ii) = 2^ii <= dob_dd;
    dob_dd = rem(dob_dd, 2^ii);
end
%Get binary no for mm
mm_bin = zeros(1,4);
for ii = 3:-1:0
    mm_bin(4 - ii) = 2^ii <= dob_mm;
    dob_mm = rem(dob_mm, 2^ii);
end
%Get binary no for yy
yy_bin = zeros(1,11);
for ii = 10:-1:0
    yy_bin(11 - ii) = 2^ii <= dob_yy;
    dob_yy = rem(dob_yy, 2^ii);
end

%Make a barcode from the nhs_id number and dob - ensure we have zeros at
%start and end
nhs_barcode = kron([0 nhs_bin dd_bin mm_bin yy_bin 0], ones(20, 4));

% Put the componenets together to build the complete form
form = [marker1 ones(50, 550);...
        ones(200, 600);
        line_scale; 
        ones(25, 600);
        line_scale;
        ones(25, 600);
        yes_no;
        ones(25, 600);
        marker2 ones(50, 500) marker4];

[r c] = size(form);

%Put the barcode/patient id box in the top write corner of the form
[bar_r bar_c] = size(nhs_barcode);
form(1:bar_r, 350:350+bar_c-1) = nhs_barcode;

%Now create a figure with the form that we can then save to a pdf

%Define size in mm of a4 paper
a4 = [210 297];

%Create figure holder
f1 = figure(...
    'Units', 'centimeters',...
    'Position', [0, 0, 21, 29.7],...
    'PaperPositionMode','auto',...
    'WindowStyle', 'normal',...
    'Color', [1, 1, 1],...
    'Visible', 'off');

%Compute bottom left coordinates (as a normalised fraction) for form
nx = c / (4*a4(1));
ny = r / (4*a4(2));

%Create axes to position form in figure
a1 = axes(...
    'Units', 'normalized',...
    'position', [(1-nx)/2, (1-nx)/2, nx, ny]);
imagesc(form); colormap(gray(256));

%Add the next to the form
text(50, 140, '\bf TamPrev Breast Density VAS Form');
view_label = {'RML', 'LML', 'RCC', 'LCC'};
app_label = {'\bf Current:', '\bf Previous:'};
for jj = 1:2
    text(15, (jj-1)*225 + 250, app_label{jj});
    for ii = 1:4
        y = 50*(ii-1) + (jj-1)*225 + 275;
        text(25, y, [view_label{ii}]);
        text(85, y, '0'); 
        text(505, y, '100');
    end
end

text(25, 735, 'Reduction of at least 10% points?' );
text(260, 735, 'Yes:' );
text(360, 735, 'No:' );

text(350, bar_r+20, ['Name: ', full_name]);%
text(350, bar_r+45, ['NHS no: ', num2str(nhs_id)]);%
text(350, bar_r+70, ['DOB: ', dob]);%
text(350, bar_r+95, ['Subject no: ', num2str(sx_id)]);%
text(350, bar_r+120, ['Current appointment: ', app_date_curr]);%
text(350, bar_r+145, ['Previous appointment: ', app_date_prev]);%
text(350, bar_r+170, 'Reader initials:          .......................');%

%Hide the axes
set(a1, 'Visible', 'off');

%Set more figure properties for printing
set(f1,...
    'PaperUnits', 'centimeters',...
    'PaperType', 'A4',...
    'PaperPosition', [0 0 21 29.7]);

%Print to pdf
form_name = [form_folder, second_name,  '_', zerostr(nhs_id, 10),'.pdf'];
print('-dpdf', '-noui', '-painters', ['-f' num2str(f1)], '-r864', form_name);

%Close the figure
close(f1);

function newstr = zerostr(num, len)
    if num >= 10^9
        newstr = num2str(num);
    else
        newstr = sprintf(['%0' num2str(len) 'u'], num);
    end
        
