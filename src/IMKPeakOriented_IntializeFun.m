function [Energy,Backbone_j_1,Flag,u0, du_i_1] = IMKPeakOriented_IntializeFun(Backbone)
%IMK PINCHING HYSTERESIS MODEL (Ibarra et al., 2005; Lignos and Krawinkler, 2011)
%##################################################################################################################################################################### 
% INPUT:--------------------------------------------------------------------------------------------------------------------------------------------------------------


% OUTPUT:--------------------------------------------------------------------------------------------------------------------------------------------------------------
% 1-  Response = Matrix of force-deformation response
%###################################################################################################################################################################### 
%% Backbone Input Parameters
%######################################################################################################################################################################


%######################################################################################################################################################################
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------	 	
Backbone_j_1.Upeak_pos  =  Backbone.Uy_pos;             
Backbone_j_1.Fpeak_pos  =  Backbone.Fy_pos;             
Backbone_j_1.Upeak_neg  = -Backbone.Uy_neg;            
Backbone_j_1.Fpeak_neg  = -Backbone.Fy_neg; 

Backbone_j_1.FLastPeak_pos  =  Backbone.Fy_pos;            
Backbone_j_1.ULastPeak_pos  =  Backbone.Uy_pos; 
Backbone_j_1.FLastPeak_neg  = -Backbone.Fy_neg;            
Backbone_j_1.ULastPeak_neg  = -Backbone.Uy_neg; 

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------	 				
Backbone_j_1.Uy_pos  =  Backbone.Uy_pos;                
Backbone_j_1.Fy_pos  =  Backbone.Fy_pos;                
Backbone_j_1.Kp_pos  =  Backbone.Kp_pos;              
Backbone_j_1.Kpc_pos  =  -Backbone.Kpc_pos;              
Backbone_j_1.Uy_neg  = -Backbone.Uy_neg;               
Backbone_j_1.Fy_neg  = -Backbone.Fy_neg;               
Backbone_j_1.Kp_neg  =  Backbone.Kp_neg;              
Backbone_j_1.Kpc_neg  =  -Backbone.Kpc_neg;              
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------	 		
Backbone_j_1.Umax_pos  =  Backbone.Umax_pos;                
Backbone_j_1.Fmax_pos  =  Backbone.Fmax_pos;
Backbone_j_1.Fres_pos  =  Backbone.Fy_pos*Backbone.FresFy_pos;
Backbone_j_1.Umax_neg  = -Backbone.Umax_neg;                
Backbone_j_1.Fmax_neg  = -Backbone.Fmax_neg;  
Backbone_j_1.Fres_neg  = -Backbone.Fy_neg*Backbone.FresFy_neg;  
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------	 		
Backbone_j_1.KrelA_pos  = Backbone.Ke;            
Backbone_j_1.KrelA_neg  = Backbone.Ke;            
Backbone_j_1.KrelB_pos  = Backbone.Ke;            
Backbone_j_1.KrelB_neg  = Backbone.Ke;  
Backbone_j_1.Krel_pos  = Backbone.Ke;  
Backbone_j_1.Krel_neg  = Backbone.Ke;  

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------		
Backbone_j_1.Kul  = Backbone.Ke;                   
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------				

Backbone_j_1.Ures_pos  = (Backbone_j_1.Fres_pos-Backbone_j_1.Fmax_pos)/Backbone_j_1.Kpc_pos+Backbone_j_1.Umax_pos;                
Backbone_j_1.Ures_neg  = (Backbone_j_1.Fres_neg-Backbone_j_1.Fmax_neg)/Backbone_j_1.Kpc_neg+Backbone_j_1.Umax_neg;               

Backbone_j_1.Ubp_pos  =  0;              
Backbone_j_1.Ubp_neg  = -0;              

Backbone_j_1.Fbp_pos  =  0;            
Backbone_j_1.Fbp_neg  = -0;  

Energy.Energy_Acc= 0.0;
Energy.Energy_Diss = 0.0; 

Flag.Failure=0;
Flag.Excursion=0;
Flag.TargetPeak=0;
Flag.Reload=0;
Flag.Unload=0;
Flag.Yield=0;
Flag.Reversal=0;

u0 =  0.0;

Backbone_j_1.Plastic_Offset_pos=0;
Backbone_j_1.Plastic_Offset_neg=0;

du_i_1=0;
%######################################################################################################################################################################
