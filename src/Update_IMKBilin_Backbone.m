function [Backbone]=Update_IMKBilin_Backbone (Kx)

global MainDirectory ProjectPath ProjectName
clc;
cd (ProjectPath) 
load(ProjectName)  
cd (MainDirectory)

if RSA_Option==2 % Maintain Fy
    
    BackboneNoPD.Ke =Kx;
    
    % Update Uy
    BackboneNoPD.Uy_pos0     = BackboneNoPD.Fy_pos0 / BackboneNoPD.Ke;
    % Update affected paramters (only Umax in this case)
    BackboneNoPD.Umax_pos0 = BackboneNoPD.Uy_pos0 + BackboneNoPD.Up_pos0;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    Kp_pos=(BackboneNoPD.Fmax_pos0-BackboneNoPD.Fy_pos0)/BackboneNoPD.Up_pos0;
    if Kp_pos>Kx
        Kp_pos=Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_pos0 = BackboneNoPD.Fy_pos0 + Kp_pos * BackboneNoPD.Up_pos0;
        BackboneNoPD.FmaxFy_pos0 = BackboneNoPD.Fmax_pos0 / BackboneNoPD.Fy_pos0;
        % Re-calulcaute Kpc and Ures based on the re-calculated Fmax
        BackboneNoPD.Kpc_pos0  = BackboneNoPD.Fmax_pos0 / BackboneNoPD.Upc_pos0;
        BackboneNoPD.Ures_pos0 = BackboneNoPD.Uy_pos0 + BackboneNoPD.Up_pos0 + (BackboneNoPD.Fmax_pos0-BackboneNoPD.Fres_pos0) / BackboneNoPD.Kpc_pos0;  
    end 

    % Update Uy
    BackboneNoPD.Uy_neg0     = BackboneNoPD.Fy_neg0 / BackboneNoPD.Ke;
    % Update affected paramters (only Umax in this case)    
    BackboneNoPD.Umax_neg0 = BackboneNoPD.Uy_neg0 + BackboneNoPD.Up_neg0;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    Kp_neg=(BackboneNoPD.Fmax_neg0-BackboneNoPD.Fy_neg0)/(BackboneNoPD.Umax_neg0-BackboneNoPD.Uy_neg0);
    if Kp_neg>Kx
        Kp_neg=Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_neg0 = BackboneNoPD.Fy_neg0 + Kp_neg * BackboneNoPD.Up_neg0;
        BackboneNoPD.FmaxFy_neg0 = BackboneNoPD.Fmax_neg0 / BackboneNoPD.Fy_neg0;
        % Re-calulcaute Kpc and Ures based on the re-calculated Fmax
        BackboneNoPD.Kpc_neg0  = BackboneNoPD.Fmax_neg0 / (BackboneNoPD.Upc_neg0);
        BackboneNoPD.Ures_neg0 = BackboneNoPD.Uy_neg0 + BackboneNoPD.Up_neg0 + (BackboneNoPD.Fmax_neg0-BackboneNoPD.Fres_neg0) / BackboneNoPD.Kpc_neg0;  
    end

    % Deduce Backbone with P-Delta 
    % Modify the strength paramters of the backbone curve to account for P-Delta
    Backbone.Fy_pos0     = BackboneNoPD.Fy_pos0 - (P*BackboneNoPD.Uy_pos0/H);
    Backbone.Fmax_pos0   = BackboneNoPD.Fmax_pos0 - P*(BackboneNoPD.Uy_pos0+BackboneNoPD.Up_pos0)/H;
    Backbone.FmaxFy_pos0 = Backbone.FmaxFy_pos0;
    Backbone.Fres_pos0   = BackboneNoPD.Fres_pos0 - (P*(BackboneNoPD.Uy_pos0+BackboneNoPD.Up_pos0+(BackboneNoPD.Fmax_pos0-BackboneNoPD.Fres_pos0)/BackboneNoPD.Kpc_pos0)/H);
    Backbone.FresFy_pos0 = Backbone.Fres_pos0 / Backbone.Fy_pos0;
    if Backbone.FresFy_pos0<0.0; Backbone.FresFy_pos0 = 0; end

    Backbone.Fy_neg0     = BackboneNoPD.Fy_neg0 - (P*BackboneNoPD.Uy_pos0/H);
    Backbone.Fmax_neg0   = BackboneNoPD.Fmax_neg0 - P*(BackboneNoPD.Uy_neg0+BackboneNoPD.Up_neg0)/H;
    Backbone.FmaxFy_neg0 = Backbone.Fmax_neg0 / Backbone.Fy_neg0;
    Backbone.Fres_neg0   = BackboneNoPD.Fres_neg0 - (P*(BackboneNoPD.Uy_neg0+BackboneNoPD.Up_neg0+(BackboneNoPD.Fmax_neg0-BackboneNoPD.Fres_neg0)/BackboneNoPD.Kpc_neg0)/H);
    Backbone.FresFy_neg0 = Backbone.Fres_neg0 / Backbone.Fy_neg0;
    if Backbone.FresFy_neg0<0.0; Backbone.FresFy_neg0 = 0; end

    Backbone.Ke =Backbone.Fy_pos0/BackboneNoPD.Uy_pos0;

    % The rest of the parameters remain the same
    Backbone.Uy_pos0   = BackboneNoPD.Uy_pos0;
    Backbone.Up_pos0   = BackboneNoPD.Up_pos0;
    Backbone.Upc_pos0  = BackboneNoPD.Upc_pos0;
    Backbone.Uu_pos0   = BackboneNoPD.Uu_pos0;
    Backbone.Uy_pos0   = BackboneNoPD.Uy_pos0;
    Backbone.Umax_pos0 = BackboneNoPD.Umax_pos0;
    
    Backbone.Uy_neg0   = BackboneNoPD.Uy_neg0;
    Backbone.Up_neg0   = BackboneNoPD.Up_neg0;
    Backbone.Upc_neg0  = BackboneNoPD.Upc_neg0;
    Backbone.Uu_neg0   = BackboneNoPD.Uu_neg0;
    Backbone.Uy_neg0   = BackboneNoPD.Uy_neg0;
    Backbone.Umax_neg0 = BackboneNoPD.Umax_neg0;   

    Backbone.Ures_pos0 = BackboneNoPD.Ures_pos0 - BackboneNoPD.Umax_pos0;  
    Backbone.Ures_neg0 = BackboneNoPD.Ures_neg0 - BackboneNoPD.Umax_neg0;  
    
    Backbone.Kpc_pos0 = (Backbone.Fmax_pos0-Backbone.Fres_pos0)/(Backbone.Ures_pos0);
    Backbone.Kpc_neg0 = (Backbone.Fmax_neg0-Backbone.Fres_neg0)/(Backbone.Ures_neg0);

    Backbone.Kp_pos0 = (Backbone.Fmax_pos0-Backbone.Fy_pos0)/(Backbone.Up_pos0);
    Backbone.Kp_neg0 = (Backbone.Fmax_neg0-Backbone.Fy_neg0)/(Backbone.Up_neg0);
    
    % Check if unloading slope (Kpc) is less than the hardening slope (Kp) and
    % modify Kul if this is the case
    if Backbone.Kpc_pos0<Backbone.Kp_pos0
        Backbone.Kpc_pos0=Backbone.Kp_pos0;
    end
    if Backbone.Kpc_neg0<Backbone.Kp_neg0
        Backbone.Kpc_neg0=Backbone.Kp_neg0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if RSA_Option==3 % Maintain Uy

    BackboneNoPD.Ke =Kx;
    
    % Update Fy
    BackboneNoPD.Fy_pos0 = BackboneNoPD.Uy_pos0 * BackboneNoPD.Ke;
    % Update affected paramters (Fmax and Fres in this case)
    BackboneNoPD.Fmax_pos0 = BackboneNoPD.FmaxFy_pos0 * BackboneNoPD.Fy_pos0;
    BackboneNoPD.Fres_pos0 = BackboneNoPD.FresFy_pos0 * BackboneNoPD.Fy_pos0;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    Kp_pos=(BackboneNoPD.Fmax_pos0-BackboneNoPD.Fy_pos0)/(BackboneNoPD.Up_pos0);
    if Kp_pos>Kx
        Kp_pos=Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_pos0 = BackboneNoPD.Fy_pos0 + Kp_pos * BackboneNoPD.Up_pos0;
        BackboneNoPD.FmaxFy_pos0 = BackboneNoPD.Fmax_pos0 / BackboneNoPD.Fy_pos0;
        % Re-calulcaute Kpc and Ures based on the re-calculated Fmax
        BackboneNoPD.Kpc_pos0  = BackboneNoPD.Fmax_pos0 / BackboneNoPD.Upc_pos0;
        BackboneNoPD.Ures_pos0 = BackboneNoPD.Uy_pos0 + BackboneNoPD.Up_pos0 + (BackboneNoPD.Fmax_pos0-BackboneNoPD.Fres_pos0) / BackboneNoPD.Kpc_pos0;  
    end 

    % Update Fy
    BackboneNoPD.Fy_neg0 = BackboneNoPD.Uy_neg0 * BackboneNoPD.Ke;
    % Update affected paramters (Fmax and Fres in this case)
    BackboneNoPD.Fmax_neg0 = BackboneNoPD.FmaxFy_neg0 * BackboneNoPD.Fy_neg0;
    BackboneNoPD.Fres_neg0 = BackboneNoPD.FresFy_neg0 * BackboneNoPD.Fy_neg0;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    Kp_neg=(BackboneNoPD.Fmax_neg0-BackboneNoPD.Fy_neg0)/(BackboneNoPD.Up_neg0);
    if Kp_neg>Kx
        Kp_neg=Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_neg0 = BackboneNoPD.Fy_neg0 + Kp_neg * BackboneNoPD.Up_neg0;
        BackboneNoPD.FmaxFy_neg0 = BackboneNoPD.Fmax_neg0 / BackboneNoPD.Fy_neg0;
        % Re-calulcaute Kpc and Ures based on the re-calculated Fmax
        BackboneNoPD.Kpc_neg0  = BackboneNoPD.Fmax_neg0 / (BackboneNoPD.Upc_neg0);
        BackboneNoPD.Ures_neg0 = BackboneNoPD.Uy_neg0 + BackboneNoPD.Up_neg0 + (BackboneNoPD.Fmax_neg0-BackboneNoPD.Fres_neg0) / BackboneNoPD.Kpc_neg0;  
    end

    % Deduce Backbone with P-Delta 
    % Modify the strength paramters of the backbone curve to account for P-Delta
    Backbone.Fy_pos0     = BackboneNoPD.Fy_pos0   - (P*BackboneNoPD.Uy_pos0/H);
    Backbone.Fmax_pos0   = BackboneNoPD.Fmax_pos0 - (P*BackboneNoPD.Umax_pos0/H);
    Backbone.FmaxFy_pos0 = BackboneNoPD.FmaxFy_pos0;
    Backbone.Fres_pos0   = BackboneNoPD.Fres_pos0 - (P*(BackboneNoPD.Umax_pos0+(BackboneNoPD.Fmax_pos0-BackboneNoPD.Fres_pos0)/BackboneNoPD.Kpc_pos0)/H);
    Backbone.FresFy_pos0 = Backbone.Fres_pos0 / Backbone.Fy_pos0;
    if Backbone.FresFy_pos0<0.0; Backbone.FresFy_pos0 = 0; end

    Backbone.Fy_neg0     = BackboneNoPD.Fy_neg0   - (P*BackboneNoPD.Uy_pos0/H);
    Backbone.Fmax_neg0   = BackboneNoPD.Fmax_neg0 - (P*BackboneNoPD.Umax_neg0/H);
    Backbone.FmaxFy_neg0 = Backbone.Fmax_neg0 / Backbone.Fy_neg0;
    Backbone.Fres_neg0   = BackboneNoPD.Fres_neg0 - (P*(BackboneNoPD.Umax_neg0+(BackboneNoPD.Fmax_neg0-BackboneNoPD.Fres_neg0)/BackboneNoPD.Kpc_neg0)/H);
    Backbone.FresFy_neg0 = Backbone.Fres_neg0 / Backbone.Fy_neg0;
    if Backbone.FresFy_neg0<0.0; Backbone.FresFy_neg0 = 0; end

    Backbone.Ke = Backbone.Fy_pos0 / BackboneNoPD.Uy_pos0;

    % The rest of the parameters remain the same
    Backbone.Uy_pos0   = BackboneNoPD.Uy_pos0;
    Backbone.Up_pos0   = BackboneNoPD.Up_pos0;
    Backbone.Upc_pos0  = BackboneNoPD.Upc_pos0;
    Backbone.Uu_pos0   = BackboneNoPD.Uu_pos0;
    Backbone.Uy_pos0   = BackboneNoPD.Uy_pos0;
    Backbone.Umax_pos0 = BackboneNoPD.Umax_pos0;
    
    Backbone.Uy_neg0   = BackboneNoPD.Uy_neg0;
    Backbone.Up_neg0   = BackboneNoPD.Up_neg0;
    Backbone.Upc_neg0  = BackboneNoPD.Upc_neg0;
    Backbone.Uu_neg0   = BackboneNoPD.Uu_neg0;
    Backbone.Uy_neg0   = BackboneNoPD.Uy_neg0;
    Backbone.Umax_neg0 = BackboneNoPD.Umax_neg0;   

    Backbone.Ures_pos0 = BackboneNoPD.Ures_pos0 - BackboneNoPD.Umax_pos0;  
    Backbone.Ures_neg0 = BackboneNoPD.Ures_neg0 - BackboneNoPD.Umax_neg0;  
    
    Backbone.Kpc_pos0 = (Backbone.Fmax_pos0-Backbone.Fres_pos0)/(Backbone.Ures_pos0);
    Backbone.Kpc_neg0 = (Backbone.Fmax_neg0-Backbone.Fres_neg0)/(Backbone.Ures_neg0);
    
    Backbone.Kp_pos0 = (Backbone.Fmax_pos0-Backbone.Fy_pos0)/(Backbone.Up_pos0);
    Backbone.Kp_neg0 = (Backbone.Fmax_neg0-Backbone.Fy_neg0)/(Backbone.Up_neg0);
    
    % Check if unloading slope (Kpc) is less than the hardening slope (Kp) and
    % modify Kul if this is the case
    if Backbone.Kpc_pos0<Backbone.Kp_pos0
        Backbone.Kpc_pos0=Backbone.Kp_pos0;
    end
    if Backbone.Kpc_neg0<Backbone.Kp_neg0
        Backbone.Kpc_neg0=Backbone.Kp_neg0;
    end
    
end