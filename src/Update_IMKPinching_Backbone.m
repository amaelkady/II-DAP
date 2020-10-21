function [Backbone]=Update_IMKPinching_Backbone (Kx)
global MainDirectory ProjectPath ProjectName
clc;
cd (ProjectPath) 
load(ProjectName)  
cd (MainDirectory)

if RSA_Option==2 % Maintain Fy
    
    BackboneNoPD.Ke = Kx;
    
    % Update Uy
    BackboneNoPD.Uy_pos  = BackboneNoPD.Fy_pos / BackboneNoPD.Ke;
    % Update affected paramters (only Umax in this case)
    BackboneNoPD.Umax_pos = BackboneNoPD.Uy_pos + BackboneNoPD.Up_pos;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    if BackboneNoPD.Kp_pos > Kx
        BackboneNoPD.Kp_pos = Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_pos   = BackboneNoPD.Fy_pos + BackboneNoPD.Kp_pos * BackboneNoPD.Up_pos;
        BackboneNoPD.FmaxFy_pos = BackboneNoPD.Fmax_pos / BackboneNoPD.Fy_pos;
        % Re-calulcaute Kpc and Ures based on the re-calculated Fmax
        BackboneNoPD.Kpc_pos  = BackboneNoPD.Fmax_pos/(BackboneNoPD.Upc_pos);
        BackboneNoPD.Ures_pos = BackboneNoPD.Uy_pos + BackboneNoPD.Up_pos + (BackboneNoPD.Fmax_pos-BackboneNoPD.Fres_pos)/BackboneNoPD.Kpc_pos; 
    end
    
    % Update Uy
    BackboneNoPD.Uy_neg  = BackboneNoPD.Fy_neg / BackboneNoPD.Ke;
    % Update affected paramters (only Umax in this case)
    BackboneNoPD.Umax_neg = BackboneNoPD.Uy_neg + BackboneNoPD.Up_neg;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    if BackboneNoPD.Kp_neg > Kx
        BackboneNoPD.Kp_neg = Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_neg   = BackboneNoPD.Fy_neg + BackboneNoPD.Kp_neg * BackboneNoPD.Up_neg;
        BackboneNoPD.FmaxFy_neg = BackboneNoPD.Fmax_neg / BackboneNoPD.Fy_neg;
        % Re-calulcaute Kpc and Ures based on the re-calculated Fmax
        BackboneNoPD.Kpc_neg  = BackboneNoPD.Fmax_neg / (BackboneNoPD.Upc_neg);
        BackboneNoPD.Ures_neg = BackboneNoPD.Uy_neg + BackboneNoPD.Up_neg + (BackboneNoPD.Fmax_neg-BackboneNoPD.Fres_neg)/BackboneNoPD.Kpc_neg;  
    end     

    % Deduce Backbone with P-Delta 
    % Modify the strength paramters of the backbone curve to account for P-Delta
    Backbone.Fy_pos     = BackboneNoPD.Fy_pos - (P*BackboneNoPD.Uy_pos/H);
    Backbone.Fmax_pos   = BackboneNoPD.Fmax_pos - P*(BackboneNoPD.Uy_pos+BackboneNoPD.Up_pos)/H;
    Backbone.FmaxFy_pos = BackboneNoPD.FmaxFy_pos;
    Backbone.Fres_pos   = BackboneNoPD.Fres_pos - (P*(BackboneNoPD.Uy_pos+BackboneNoPD.Up_pos+(BackboneNoPD.Fmax_pos-BackboneNoPD.Fres_pos)/BackboneNoPD.Kpc_pos)/H);
    Backbone.FresFy_pos = Backbone.Fres_pos / Backbone.Fy_pos;
    if Backbone.FresFy_pos<0.0; Backbone.FresFy_pos=0; end

    Backbone.Fy_neg     = BackboneNoPD.Fy_neg - (P*BackboneNoPD.Uy_pos/H);
    Backbone.Fmax_neg   = BackboneNoPD.Fmax_neg - P*(BackboneNoPD.Uy_neg+BackboneNoPD.Up_neg)/H;
    Backbone.FmaxFy_neg = Backbone.Fmax_neg / Backbone.Fy_neg;
    Backbone.Fres_neg   = BackboneNoPD.Fres_neg - (P*(BackboneNoPD.Uy_neg+BackboneNoPD.Up_neg+(BackboneNoPD.Fmax_neg-BackboneNoPD.Fres_neg)/BackboneNoPD.Kpc_neg)/H);
    Backbone.FresFy_neg = Backbone.Fres_neg / Backbone.Fy_neg;
    if Backbone.FresFy_neg<0.0; Backbone.FresFy_neg=0; end

    Backbone.Ke =Backbone.Fy_pos/BackboneNoPD.Uy_pos;

    % The rest of the parameters remain the same
    Backbone.Up_pos   = BackboneNoPD.Up_pos;
    Backbone.Upc_pos  = BackboneNoPD.Upc_pos;
    Backbone.Uu_pos   = BackboneNoPD.Uu_pos;
    Backbone.Uy_pos   = BackboneNoPD.Uy_pos;
    Backbone.Umax_pos = BackboneNoPD.Umax_pos;

    Backbone.Up_neg   = BackboneNoPD.Up_neg;
    Backbone.Upc_neg  = BackboneNoPD.Upc_neg;
    Backbone.Uu_neg   = BackboneNoPD.Uu_neg;
    Backbone.Uy_neg   = BackboneNoPD.Uy_neg;
    Backbone.Umax_neg = BackboneNoPD.Umax_neg;

    Backbone.Kp_pos = (Backbone.FmaxFy_pos-1) * Backbone.Fy_pos/(Backbone.Up_pos);
    Backbone.Kp_neg = (Backbone.FmaxFy_neg-1) * Backbone.Fy_neg/(Backbone.Up_neg);

    Backbone.Ures_pos = BackboneNoPD.Ures_pos - BackboneNoPD.Umax_pos;  
    Backbone.Ures_neg = BackboneNoPD.Ures_neg - BackboneNoPD.Umax_neg;   

    Backbone.Kpc_pos = (Backbone.Fmax_pos-Backbone.Fres_pos)/(Backbone.Ures_pos);
    Backbone.Kpc_neg = (Backbone.Fmax_neg-Backbone.Fres_neg)/(Backbone.Ures_neg);

    Backbone.Kp_pos = (Backbone.Fmax_pos-Backbone.Fy_pos)/(Backbone.Up_pos);
    Backbone.Kp_neg = (Backbone.Fmax_neg-Backbone.Fy_neg)/(Backbone.Up_neg);
    
    % Check if unloading slope (Kpc) is less than the hardening slope (Kp) and
    % modify Kul if this is the case
    if Backbone.Kpc_pos<Backbone.Kp_pos
        Backbone.Kpc_pos=Backbone.Kp_pos;
    end
    if Backbone.Kpc_neg<Backbone.Kp_neg
        Backbone.Kpc_neg=Backbone.Kp_neg;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if RSA_Option==3 % Maintain Uy

    BackboneNoPD.Ke =Kx;
    
    % Update Fy
    BackboneNoPD.Fy_pos  = BackboneNoPD.Uy_pos * BackboneNoPD.Ke;
    % Update affected paramters (Fmax and Fres in this case)
    BackboneNoPD.Fmax_pos = BackboneNoPD.FmaxFy_pos * BackboneNoPD.Fy_pos;
    BackboneNoPD.Fres_pos = BackboneNoPD.FresFy_pos * BackboneNoPD.Fy_pos;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    BackboneNoPD.Kp_pos = (BackboneNoPD.Fmax_pos-BackboneNoPD.Fy_pos)/(BackboneNoPD.Umax_pos-BackboneNoPD.Uy_pos);
    if BackboneNoPD.Kp_pos > Kx
        BackboneNoPD.Kp_pos = Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_pos   = BackboneNoPD.Fy_pos + BackboneNoPD.Kp_pos * BackboneNoPD.Up_pos;
        BackboneNoPD.FmaxFy_pos = BackboneNoPD.Fmax_pos / BackboneNoPD.Fy_pos;
        % Re-calulcaute Kpc and Ures based on the re-calculated Fmax
        BackboneNoPD.Kpc_pos  = BackboneNoPD.Fmax_pos/(BackboneNoPD.Upc_pos);
        BackboneNoPD.Ures_pos = BackboneNoPD.Uy_pos+BackboneNoPD.Up_pos+(BackboneNoPD.Fmax_pos-BackboneNoPD.Fres_pos)/BackboneNoPD.Kpc_pos; 
    end 

    % Update Fy    
    BackboneNoPD.Fy_neg  = BackboneNoPD.Uy_neg * BackboneNoPD.Ke;
    % Update affected paramters (Fmax and Fres in this case)
    BackboneNoPD.Fmax_neg = BackboneNoPD.FmaxFy_neg * BackboneNoPD.Fy_neg;
    BackboneNoPD.Fres_neg = BackboneNoPD.FresFy_neg * BackboneNoPD.Fy_neg;
    % Check if hardening slope (Kp) exceeds the elastic slope (Ke) and
    % modify Kp if this is the case
    BackboneNoPD.Kp_neg = (BackboneNoPD.Fmax_neg-BackboneNoPD.Fy_neg)/(BackboneNoPD.Umax_neg-BackboneNoPD.Uy_neg);
    if BackboneNoPD.Kp_neg > Kx
        BackboneNoPD.Kp_neg = Kx;
        % Re-calulcaute Fmax using the updated Kp and while maintaining theta_p 
        BackboneNoPD.Fmax_neg   = BackboneNoPD.Fy_neg + BackboneNoPD.Kp_neg * BackboneNoPD.Up_neg;
        BackboneNoPD.FmaxFy_neg = BackboneNoPD.Fmax_neg / BackboneNoPD.Fy_neg;
        BackboneNoPD.Kpc_neg  = BackboneNoPD.Fmax_neg / BackboneNoPD.Upc_neg;
        BackboneNoPD.Ures_neg = BackboneNoPD.Uy_neg + BackboneNoPD.Up_neg + (BackboneNoPD.Fmax_neg-BackboneNoPD.Fres_neg)/BackboneNoPD.Kpc_neg;  
    end     

    % Deduce Backbone with P-Delta 
    % Modify the strength paramters of the backbone curve to account for P-Delta
    Backbone.Fy_pos     = BackboneNoPD.Fy_pos-(P*BackboneNoPD.Uy_pos/H);
    Backbone.Fmax_pos   = BackboneNoPD.Fmax_pos-P*(BackboneNoPD.Uy_pos+BackboneNoPD.Up_pos)/H;
    Backbone.FmaxFy_pos = BackboneNoPD.FmaxFy_pos;
    Backbone.Fres_pos   = BackboneNoPD.Fres_pos-(P*(BackboneNoPD.Uy_pos+BackboneNoPD.Up_pos+(BackboneNoPD.Fmax_pos-BackboneNoPD.Fres_pos)/BackboneNoPD.Kpc_pos)/H);
    Backbone.FresFy_pos = Backbone.Fres_pos/Backbone.Fy_pos;
    if Backbone.FresFy_pos<0.0; Backbone.FresFy_pos=0; end

    Backbone.Fy_neg     = BackboneNoPD.Fy_neg - (P*BackboneNoPD.Uy_pos/H);
    Backbone.Fmax_neg   = BackboneNoPD.Fmax_neg - P*(BackboneNoPD.Uy_neg+BackboneNoPD.Up_neg)/H;
    Backbone.FmaxFy_neg = Backbone.Fmax_neg / Backbone.Fy_neg;
    Backbone.Fres_neg   = BackboneNoPD.Fres_neg - (P*(BackboneNoPD.Uy_neg+BackboneNoPD.Up_neg+(BackboneNoPD.Fmax_neg-BackboneNoPD.Fres_neg)/BackboneNoPD.Kpc_neg)/H);
    Backbone.FresFy_neg = Backbone.Fres_neg / Backbone.Fy_neg;
    if Backbone.FresFy_neg<0.0; Backbone.FresFy_neg=0; end

    Backbone.Ke = Backbone.Fy_pos / BackboneNoPD.Uy_pos;

    % The rest of the parameters remain the same
    Backbone.Up_pos   = BackboneNoPD.Up_pos;
    Backbone.Upc_pos  = BackboneNoPD.Upc_pos;
    Backbone.Uu_pos   = BackboneNoPD.Uu_pos;
    Backbone.Uy_pos   = BackboneNoPD.Uy_pos;
    Backbone.Umax_pos = BackboneNoPD.Umax_pos;

    Backbone.Up_neg   = BackboneNoPD.Up_neg;
    Backbone.Upc_neg  = BackboneNoPD.Upc_neg;
    Backbone.Uu_neg   = BackboneNoPD.Uu_neg;
    Backbone.Uy_neg   = BackboneNoPD.Uy_neg;
    Backbone.Umax_neg = BackboneNoPD.Umax_neg;

    Backbone.Kp_pos = (Backbone.FmaxFy_pos-1)*Backbone.Fy_pos/(Backbone.Up_pos);
    Backbone.Kp_neg = (Backbone.FmaxFy_neg-1)*Backbone.Fy_neg/(Backbone.Up_neg);

    Backbone.Ures_pos = BackboneNoPD.Ures_pos - BackboneNoPD.Umax_pos;  
    Backbone.Ures_neg = BackboneNoPD.Ures_neg - BackboneNoPD.Umax_neg;   

    Backbone.Kpc_pos = (Backbone.Fmax_pos-Backbone.Fres_pos)/(Backbone.Ures_pos);
    Backbone.Kpc_neg = (Backbone.Fmax_neg-Backbone.Fres_neg)/(Backbone.Ures_neg);
    
    Backbone.Kp_pos = (Backbone.Fmax_pos-Backbone.Fy_pos)/(Backbone.Up_pos);
    Backbone.Kp_neg = (Backbone.Fmax_neg-Backbone.Fy_neg)/(Backbone.Up_neg);
    
    % Check if unloading slope (Kpc) is less than the hardening slope (Kp) and
    % modify Kul if this is the case
    if Backbone.Kpc_pos0<Backbone.Kp_pos0
        Backbone.Kpc_pos0=Backbone.Kp_pos0;
    end
    if Backbone.Kpc_neg0<Backbone.Kp_neg0
        Backbone.Kpc_neg0=Backbone.Kp_neg0;
    end
end
