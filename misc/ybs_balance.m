%%
daily_rate = 1.0529^(1/365) - 1;
balance = 139930.00;
offset = 0;
tot_interest = 0;
figure;
change_balance = zeros(85, 1);
change_offset = zeros(85, 1);

%%%
 change_balance(16) = -455.19;
 change_balance(23) = -773.11;
 change_balance(51) = -773.11;
 change_balance(83) = -773.11;
%%%
 change_offset(9) = 837;
 change_offset(35) = 50;
 change_offset(37) = 15000;
 change_offset(48) = 1000;
 change_offset(49) = 1000;
 change_offset(51) = 1000;
 change_offset(52) = 1000;
 change_offset(56) = 2000;
 change_offset(62) = 450;
 change_offset(81) = 2000;
%%%

for i = 1:85
    balance = balance + change_balance(i);
    offset = offset + change_offset(i);
    interest = (balance - offset) * daily_rate;
    tot_interest = tot_interest + interest;
    balance = balance + interest;
    hold on;
    plot(i, balance, 'bx');
    %plot(i, offset, 'rx');
end