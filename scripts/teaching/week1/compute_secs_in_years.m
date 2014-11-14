function total_secs = compute_secs_in_years(num_years)
%Function to compute the total number of seconds in any given number of
%years
num_days = 365;
num_hours = 24;
num_mins = 60;
num_secs = 60;
total_secs = num_years*num_days*num_hours*num_mins*num_secs;