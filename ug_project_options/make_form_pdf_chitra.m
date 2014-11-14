function form_name = make_form_pdf_chitra(form_type, form_idx, reader_id, form_folder, debug_flag)

%Notes - as of 4th November 2009, marker 4 is used in place of marker 3 to
%distinguish between copies of the form in which the order of lines has
%changed

if nargin < 5
    debug_flag = 0;
end

%load in markers (marker1, marker2, marker3)
load components.mat marker*

% Build lines
line_scale = true(100, 600);
line_scale(49:51, 101:500) = 0;

%Work out form type from code
switch form_type
    case 1
        form_type_bin = [0 0 1];
        im_code = 'I';
        experiment_type = 'Reading Model I: Individual Assessment';
        
    case 2
        form_type_bin = [0 1 0];
        im_code = 'P';
        experiment_type = 'Reading Model II: Paired Assessment';
        
end
%Get image numbers
start_im = 5*(form_idx-1);

%Get binary no for reader_id
reader_copy = reader_id;
reader_bin = zeros(1,5);
for ii = 4:-1:0
    reader_bin(5 - ii) = 2^ii <= reader_copy;
    reader_copy = rem(reader_copy, 2^ii);
end

%Get binary no for form_id
form_copy = form_idx;
form_bin = zeros(1,5);
for ii = 4:-1:0
    form_bin(5 - ii) = 2^ii <= form_copy;
    form_copy = rem(form_copy, 2^ii);
end

%Make a barcode from the nhs_id number and dob - ensure we have zeros at
%start and end
form_barcode = kron([0 reader_bin form_type_bin form_bin 0], ones(30, 4));

% Put the componenets together to build the complete form
form = [marker1 ones(50, 550);...
        ones(200, 600);
        line_scale; line_scale; line_scale; line_scale; line_scale;
        ones(50, 600);
        marker2 ones(50, 500) marker4];

[r c] = size(form);

%Put the barcode/patient id box in the top write corner of the form
[bar_r bar_c] = size(form_barcode);
form(1:bar_r, 350:350+bar_c-1) = form_barcode;

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
text(50, 140, 'Breast Density VAS Form');
for ii = 1:5
    text(45, 200 + ii*100, [im_code num2str(start_im + ii)]);
    text(85, 200 + ii*100, '0'); 
    text(505, 200 + ii*100, '100');
end

text(350, bar_r+20, ['Reader ID: ', num2str(reader_id)]);%
text(350, bar_r+45, experiment_type);%
text(350, bar_r+70, ['Images: ', num2str(start_im+1) ' to ' num2str(start_im+5)]);%
text(350, bar_r+95, 'Reader initials:   .......................');%

%Hide the axes
set(a1, 'Visible', 'off');

%Set more figure properties for printing
set(f1,...
    'PaperUnits', 'centimeters',...
    'PaperType', 'A4',...
    'PaperPosition', [0 0 21 29.7]);

if debug_flag
    set(f1, 'visible', 'on');
else
    %Print to pdf
    form_name = [form_folder 'reader' num2str(reader_id)  '_' num2str(form_type)  '_' zerostr(form_idx, 2) '.pdf'];
    print('-dpdf', '-noui', '-painters', ['-f' num2str(f1)], '-r864', form_name);

    %Close the figure
    close(f1);
end
