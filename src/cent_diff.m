%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SDOF solution using Central Difference Method
%
% Developed by: Dimitrios G. Lignos, Ph.D.
% Date: 08/21/2011
% Last Modified: 08/21/2011
%
% Input Parameters:
% mass =  1.0
% Period T 
% deltat = time step (equal to time step of the earthquake
% xi = damping ratio
% GMotion
% g: multiples of g in proper units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Amax] =  cent_diff(p,period, deltat, xi, g)

 mass = 1.0;
 
%% Read the ground motion file and put everything in a column vector

p = p*g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have an array p() which contains L acceleration values 
% at an interval of dt sec
 
 % omega = natural frequency
 % stf = stiffness
 % damping = "c"
   
    omega = 2.0*pi/period;
    stf   = mass*omega*omega;
    damp = 2.0*mass*omega*xi;


 %% Initial Conditions of Central Difference Method
    a_prev = p(1)/mass;
    u_prev = 0.0;
    v_prev = 0.0;
    u_pp = a_prev*deltat*deltat/2.0; 
    
    consA = mass/deltat/deltat - damp/(2*deltat);
    consB = stf - 2*mass/(deltat*deltat);
    kbar = mass/deltat/deltat + damp/(2*deltat);
    time(1) = 0.0;

 % u_curr = displacement at step (i+1)
 % v_curr = velocity at step (i+1)
 % a_curr = acceleration at step (i+1)
%% Calculations for time step i 
    for i = 2:length(p)
      p_curr = p(i) - consA*u_pp - consB*u_prev;
      u_curr = p_curr/kbar;
      v_curr = (u_curr - u_pp)/(2*deltat);
      a_curr = (u_curr - 2*u_prev + u_pp)/(deltat^2);
    
      dis(i) = 1.*u_curr;
% Spectral Ordinate
      aa(i) = ((2*pi/period)^2)*u_curr;
      
      vv(i) = (2*pi/period)*u_curr;
    % store current values into (i-1) for next step
    
      u_pp = u_prev;
      u_prev = u_curr;
      v_prev = v_curr;
      a_prev = a_curr;
      
      time(i) = time(i-1) + deltat;
  end
        
    Amax = max(abs(aa));