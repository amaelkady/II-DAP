function [uii,vii,aii] = NewmarkSingleIncrNLX(ki,C,M,gamma,beta,dt,ui,vi,fi,agi,agii)
%######################################################################################################################################################################
% SDOF Linear Dynamic System Newmark Integrator Single Timestep 
%#####################################################################################################################################################################
% INPUT:--------------------------------------------------------------------------------------------------------------------------------------------------------------
% 1  - ki = stiffness of SDOF at time step i
% 2  - C = damping coefficient of SDOF, C = 2 * m * zeta * omega_n
% 3  - M = mass of SDOF
% 4  - gamma = Newmark parameters = 1/2 for Linear Acceleration and Average Acceleration (gamma >= 1/2)
% 5  - beta = Newmark parameters = 1/6 for Linear Acceleration Method and 1/4 for Average Acceleration Method (beta >= 0)
% 6  - dt = time increment for current time step
% 7  - ui = relative displacment of SDOF at time step i
% 8  - vi = relative velocity of SDOF at time step i
% 9  - fi = force on SDOF at time step i, obtained from the system type
% 10 - agi = ground acceleration at time step i
% 11 - agii = ground acceleration at time step i+1
% OUTPUT:--------------------------------------------------------------------------------------------------------------------------------------------------------------
% 1 - uii = relative displacement of SDOF at time step i+1
% 2 - vii = relative velocity of SDOF at time step i+1
% 3 - aii = relative acceleration of SDOF at time step i+1
%######################################################################################################################################################################
% see Chopra's Dynamics of Structures page 177
%######################################################################################################################################################################
if dt==0
    
    uii = 0.0;
    vii = 0.0;
    aii = 0.0;   
    
else

    pi =  -M * agi;
    pii = -M * agii;
    dp = pii-pi;
        
    ai = (pi-C*vi-fi)/M;

    dk_eff = ki + C *  gamma/dt/beta + M /dt^2/beta ;
    dp_eff = dp + M * (vi/beta/dt + ai/2/beta) + C * (vi*gamma/beta + ai*dt*(gamma/(2*beta)-1));

    du=dp_eff/dk_eff;
    dv=du*gamma/beta/dt - vi*gamma/beta + ai *dt*(1-gamma/2/beta);
    da=du/beta/dt^2 - vi/beta/dt - ai/2/beta;

    uii = ui + du;
    vii = vi + dv;
    aii = ai + da;
        
end
%######################################################################################################################################################################
end

