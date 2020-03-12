function F = compute_f_groupSixty(t,Frmax,Fymax,amiapredator,pr,vr,Er,py,vy,Ey)


% PLEASE FILL OUT THE INFORMATION BELOW WHEN YOU SUBMIT YOUR CODE
%% Test time and place: Enter the time and room for your test here 
% Group members: Benjamin Joshua Mason Serdar


%   t: Time
%   Frmax: Max force that can act on the predator
%   Fymax: Max force that can act on the prey
%   amiapredator: Logical variable - if amiapredator is true,
%   the function must compute forces acting on a predator.
%   If false, code must compute forces acting on a prey.
%   pr - 2D vector with current position of predator eg pr = [x_r;y_r]
%   vr - 2D vector with current velocity of predator eg vr= [vx_r;vy_r]
%   Er - energy remaining for predator
%   py - 2D vector with current position of prey py = [x_prey;y_prey]
%   vy - 2D vector with current velocity of prey py = [vx_prey;vy_prey]
%   Ey - energy remaining for prey
%   F - 2D vector specifying the force to be applied to the object
%   that you wish to control F = [Fx;Fy]
%   The direction of the force is arbitrary, but if the
%   magnitude you specify exceeds the maximum allowable
%   value its magnitude will be reduced to this value
%   (without changing direction)

    g = 9.81;
    mr = 100; % Mass of predator, in kg
    my = 10.; % Mass of prey, in kg
    predator_crash_limit = 15; % Predator max landing speed to survive
    prey_crash_limit = 8; % Prey max landing speed to survive
    Max_fuel_r = 500000; % Max stored energy for predator
    Max_fuel_y = 50000;  % Max stored energy for prey

  if (amiapredator)
%###% Code to compute the force to be applied to the predator#############
    velI=min(.05*norm(py-pr),5); %how far ahead to look
    Ftp=((py+vy*velI)-(pr+vr*velI*.9)); %go towards future positions
    Ftp=.5*Frmax*Ftp/norm(Ftp);
    if t<=5
        Flaunch=[0;1*g*mr];
    else
        Flaunch=[0;0];
    end
    if(Er>Max_fuel_r/5 || t>200) %if don't need to refuel
            F=Ftp+Flaunch+[0;g*mr];  %add all forces together
    else
        F=[0;0]; %if high up just fall
        if(pr(2)<500) %low enough to start stopping
            F=Frmax*-vr*.07; %extra drag force slows to safe falling rate
        end
    end
  else %prey, not a predator
%###% Code to compute the force to be applied to the prey################# 
    Fgroundy=[0;160*Fymax/(py(2)+.02)^5];
        theta=70;
        dir=[cos(-theta*pi/180)*vr(1)-sin(-theta*pi/180)*vr(2);sin(-theta*pi/180)*vr(1)+cos(-theta*pi/180)*vr(2)];
        Fescy=Fymax*dir/(norm(vr)+.0001)*(.01+.6-min(norm(pr-py)^2/100000,.6));
        Fescy=Fescy/(norm(Fescy)+.001)*min(norm(Fescy),Fymax*.9);
    F=Fgroundy+Fescy+[0;my*g];
  end %end prey, not a predator
end