%Script showing examples of how data can be stored in forms other than
%standard arrays
%% Cells - creating
%Useful when we want to store data of different sizes, lengths etc together

%The names in my offices
names_c = {'Michael', 'Andy', 'Phil'}; %Note the use of curly braces {, instead of square brackets [
names_r = {'Michael'; 'Andy'; 'Phil'};
%names_c is 1x3 cell array, names_r is a 3x1 array - note how ; and , act
%the same as in standard arrays

%We can also pre-allocate cells
new_cell = cell(5, 3); %A 5x3 cell array, each cell is initialised to an empty array
new_cell{1, 2} = 'Hello';
new_cell{3, 1} = 1;
new_cell{7} = [1 2 3 4];
%We can index into the array just like regular arrays, and if we want,
%store different data types in different cells (although I wouldn't usually
%use cell in this way)

%Note that {} returns the data contained in each cell. We can also use ()
%to return the cells themselves
a_string = new_cell{1, 2}; %Get the data from the cell in the 1st row, 2nd column - a_string is 1x5 char array
a_cell = new_cell(1, 2); %Gets the cell in the 1st row, 2nd column - a_cell is a 1x1 cell (which happens to contain a 1x5 string)

%We can use () to swap parts of one cell array to another - in exactly the
%same way swapped numerical array elements in the first week

%% Structures

%The function dir is a good example that uses a structure array to collate
%information on a set of files
file_list = dir();

%We can index these like regular arrays, this returns a single element of
%the structure which will have a set of fields. Every element will have the
%same fields, even if they are empty
file_list(3,1);
file_list(:);

%We can access the fields within each element using .
f3 = file_list(3).name;

%We can also create structures explicity like this
my_struc(1).firstname = 'Mike';
my_struc(1).surname = 'Berks';
my_struc(1).age = 30;

my_struc(2).firstname = 'Phil';
my_struc(2).surname = 'Tresadern';
my_struc(2).age = 33;

%Or like this
my_struct2 = struct('firstname', {'Mike'; 'Phil'}, 'surname', {'Berks'; 'Tresadern'}, 'age', [30; 33]); %Note the use of cells, describing the values to fill each cell with

%We'll return to structures later in the course








