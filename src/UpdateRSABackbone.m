function [Backbone] = UpdateRSABackbone (TT,P_original,H_original,SC_original)
global MainDirectory ProjectPath ProjectName

cd (ProjectPath)
load(ProjectName) ;
cd (MainDirectory)

Kx=M/(TT/2/pi)^2;
Stability_Coeff=P_original/H_original/Kx;
H=P_original/Stability_Coeff/Kx;

%%
if SystemModelstatus==1
    delta_PDeltaFail_pos = 0;
    delta_PDeltaFail_neg = 0;
    
    % Modify the paramters of the backbone curve (without P-Delta) based on
    % the new Ke
    BackboneNoPD.Ke = Kx;
    
    % Deduced the bilinear curve to account for P-Delta
    BackbonePD.Ke=BackboneNoPD.Ke-P/H;

    Backbone=BackbonePD;
    
    % Get Failure displacement
    if Backbone.Ke < 0.0
        delta_PDeltaFail = 0.0;
    else
        delta_PDeltaFail = 999;
    end
end

%%
if SystemModel_Option==2
    delta_PDeltaFail_pos=0;
    delta_PDeltaFail_neg=0;
    
    % Modify the paramters of the backbone curve (without P-Delta) based on
    % the new Ke
    BackboneNoPD.Ke = Kx;
    if RSA_Option==2
        BackboneNoPD.dy=BackboneNoPD.Fy / BackboneNoPD.Ke; % maintain original Fy and deduced new dy
    elseif RSA_Option==3
        BackboneNoPD.Fy=BackboneNoPD.dy * BackboneNoPD.Ke; % maintain original dy and deduced new Fy        
    end    
    BackboneNoPD.Kp = BackboneNoPD.Kp; % maintain original Kp
    if BackboneNoPD.Kp>BackboneNoPD.Ke % limit Kp to Ke in case it exceeds it
        BackboneNoPD.Kp=BackboneNoPD.Ke;
    end
    
    % Deduced the bilinear curve to account for P-Delta
    BackbonePD.dy=BackboneNoPD.dy;
    BackbonePD.Ke=BackboneNoPD.Ke-P/H;
    BackbonePD.Kp=BackboneNoPD.Kp-P/H;
    BackbonePD.Fy=BackboneNoPD.Fy-P*BackboneNoPD.dy/H;

    Backbone=BackbonePD;
    
    % Get Failure displacement
    if Backbone.Ke < 0.0
        delta_PDeltaFail=0.0;
    else
        if Backbone.Kp < 0.0
            delta_PDeltaFail = Backbone.dy + Backbone.Fy / abs(Backbone.Kp);
        else
            delta_PDeltaFail=999;
        end
    end
end

%%
if SystemModel_Option==3
    delta_PDeltaFail_pos=0;
    delta_PDeltaFail_neg=0;

    [Backbone]=Update_IMKBilin_Backbone (Kx);
    
    if Backbone.Fres_pos0 < 0.0
        delta_PDeltaFail_pos = Backbone.Uy_pos0 + Backbone.Up_pos0 + Backbone.Fmax_pos0/Backbone.Kpc_pos0;
    else
        delta_PDeltaFail_pos = Backbone.Uu_pos0;  
    end
    
    if Backbone.Fres_neg0 < 0.0
        delta_PDeltaFail_neg = Backbone.Uy_neg0 + Backbone.Up_neg0 + Backbone.Fmax_neg0/Backbone.Kpc_neg0;
    else
        delta_PDeltaFail_neg= Backbone.Uu_neg0;  
    end   
end

%%
if SystemModel_Option==4 || SystemModel_Option == 5
    delta_PDeltaFail_pos=0;
    delta_PDeltaFail_neg=0;

    [Backbone]=Update_IMKPinching_Backbone (Kx);

    if BackbonePD.Fres_pos < 0.0
        delta_PDeltaFail_pos = BackbonePD.Uy_pos + BackbonePD.Up_pos + BackbonePD.Fmax_pos/BackbonePD.Kpc_pos;
    else
        delta_PDeltaFail_pos = BackbonePD.Uu_pos;  
    end
    
    if BackbonePD.Fres_neg < 0.0
        delta_PDeltaFail_neg = BackbonePD.Uy_neg + BackbonePD.Up_neg + BackbonePD.Fmax_neg/BackbonePD.Kpc_neg;
    else
        delta_PDeltaFail_neg= BackbonePD.Uu_neg;  
    end  
    
end

%%
if SystemModel_Option==6
    delta_PDeltaFail_pos=0;
    delta_PDeltaFail_neg=0;

    [Backbone]=Update_IMKSelfCentering_Backbone (Kx);
    
    if Backbone.Fres_pos0 < 0.0
        delta_PDeltaFail_pos = Backbone.Uy_pos0 + Backbone.Up_pos0 + Backbone.Fmax_pos0/Backbone.Kpc_pos0;
    else
        delta_PDeltaFail_pos = Backbone.Uu_pos0;  
    end
    
    if Backbone.Fres_neg0 < 0.0
        delta_PDeltaFail_neg = Backbone.Uy_neg0 + Backbone.Up_neg0 + Backbone.Fmax_neg0/Backbone.Kpc_neg0;
    else
        delta_PDeltaFail_neg= Backbone.Uu_neg0;  
    end   
end
%%
cd (ProjectPath)
pause(0.1)
save(ProjectName,'Backbone','delta_PDeltaFail','delta_PDeltaFail_pos','delta_PDeltaFail_neg','-append')
pause(0.5)
cd (MainDirectory)