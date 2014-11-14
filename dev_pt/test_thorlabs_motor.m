function test_thorlabs_motor()
%% Test motor control under ActiveX
clc;

h = actxcontrol('MGMOTOR.MGMotorCtrl.1');

SN = 83837992;
set(h,'HWSerialNum', SN);

% Start Control
h.StartCtrl;
h.Identify;

pause(5); % waiting for the GUI to load up;

% Associate events with corresponding handlers
h.registerevent({'MoveComplete' 'MoveCompleteHandler'});

timeout = 10; % timeout for waiting the move to be completed
%h.MoveJog(0,1); % Jog
 
% Move a absolute distance
h.SetAbsMovePos(0, 7);
h.MoveAbsolute(0, 1==0);
 
t0 = clock;
while(etime(clock, t0) < timeout) 
% wait while the motor is active; timeout to avoid dead loop
    s = h.GetStatusBits_Bits(0);
    if ~IsMoving(s)
      pause(2); % pause 2 seconds;
      h.MoveHome(0, 0);
      disp('Home Started!');
      break;
    end
end



function MoveCompleteHandler(varargin)
disp('Move complete!');


function r = IsMoving(StatusBits)
% Read StatusBits returned by GetStatusBits_Bits method and determine if
% the motor shaft is moving; Return 1 if moving, return 0 if stationary
r = bitget(abs(StatusBits),5) || bitget(abs(StatusBits),6);


