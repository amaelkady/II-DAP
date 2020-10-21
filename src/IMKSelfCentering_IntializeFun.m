%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                          %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%     SMOOTH IMK MODEL     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                          %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Written by: Ahmed ELkady %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Date Created:  25 April  2017  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Last Modified: 28 August 2017  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This Code Reads the IMK Input Paramters and Returns Ri and Mi at a given Step

function [BackbonePos_Original, BackboneNeg_Original, BackbonePos_j_1, BackboneNeg_j_1, Beta_j_1, Constants, Flags, InitialValues, Smoothparamters, Ri_1, Mi_1, Di_1, Energy_Excrsni_1,Energy_Excrsn, Energy_Rev, Energy_total]=IMKSelfCentering_IntializeFun (Backbone)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  EXTRACT BACKBONE VALUES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ke=Backbone.Ke;

Up_pos0=Backbone.Up_pos0;
Upc_pos0=Backbone.Upc_pos0;
Uu_pos0=Backbone.Uu_pos0;
Fy_pos0=Backbone.Fy_pos0;
FmaxFy_pos0=Backbone.FmaxFy_pos0;
FresFy_pos0=Backbone.FresFy_pos0;
Fll_pos0=Backbone.Fll_pos0;
Ull_pos0=Backbone.Ull_pos0;
Kul_pos0=Backbone.Kul_pos0;

Up_neg0=Backbone.Up_neg0;
Upc_neg0=Backbone.Upc_neg0;
Uu_neg0=Backbone.Uu_neg0;
Fy_neg0=Backbone.Fy_neg0;
FmaxFy_neg0=Backbone.FmaxFy_neg0;
FresFy_neg0=Backbone.FresFy_neg0;
Fll_neg0=Backbone.Fll_neg0;
Ull_neg0=Backbone.Ull_neg0;
Kul_neg0=Backbone.Kul_neg0;

LAMBDA_S=Backbone.LAMBDA_S;
LAMBDA_C=Backbone.LAMBDA_C;
LAMBDA_K=Backbone.LAMBDA_K;
c_S=Backbone.c_S;
c_C=Backbone.c_C;
c_K=Backbone.c_K;

D_pos=Backbone.D_pos;
D_neg=Backbone.D_neg;

n=Backbone.n;
Roffset=Backbone.Roffset;
LAMBDA_F=Backbone.LAMBDA_F;
c_F=Backbone.c_F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  INITIAL VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   PRE-CALCAULTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deduced parameters
Ref_Energy_S=LAMBDA_S*Fy_pos0;
Ref_Energy_C=LAMBDA_C*Fy_pos0;
Ref_Energy_K=LAMBDA_K*Fy_pos0;
Ref_Energy_F=LAMBDA_F*Fy_pos0;

Fres_pos0=Backbone.Fres_pos0;
Fres_neg0=Backbone.Fres_neg0;

       Fmax_pos0 = FmaxFy_pos0*Fy_pos0;
    Uy_pos0 = Fy_pos0/Ke;
  Umax_pos0 = Uy_pos0+Up_pos0;
    Kp_pos0 = (Fmax_pos0-Fy_pos0)/(Up_pos0);
   Kpc_pos0 = Backbone.Kpc_pos0;
 FyProject_pos0 = Fmax_pos0-Kp_pos0 *Umax_pos0;
FmaxProject_pos0 = Fmax_pos0+Kpc_pos0*Umax_pos0;
    
       Fmax_neg0 = FmaxFy_neg0*Fy_neg0;
    Uy_neg0 = Fy_neg0/Ke;
  Umax_neg0 = Uy_neg0+Up_neg0;
    Kp_neg0 = (Fmax_neg0-Fy_neg0)/(Up_neg0);
   Kpc_neg0 = Backbone.Kpc_neg0;
 FyProject_neg0 = Fmax_neg0-Kp_neg0 *Umax_neg0;
FmaxProject_neg0 = Fmax_neg0+Kpc_neg0*Umax_neg0;

            K_j_1 = Ke;
    Uy_pos_j_1 = Uy_pos0;
  Umax_pos_j_1 = Umax_pos0;
    Kp_pos_j_1 = Kp_pos0;
   Kpc_pos_j_1 = Kpc_pos0;
    Uu_pos_j_1 = Uu_pos0;
        Fy_pos_j_1 = Fy_pos0;
 FyProject_pos_j_1 = FyProject_pos0;
       Fmax_pos_j_1 = Fmax_pos0;    
FmaxProject_pos_j_1 = FmaxProject_pos0;

    Uy_neg_j_1 = Uy_neg0;
  Umax_neg_j_1 = Umax_neg0;
    Uu_neg_j_1 = Uu_neg0;
    Kp_neg_j_1 = Kp_neg0;
   Kpc_neg_j_1 = Kpc_neg0;
        Fy_neg_j_1 = Fy_neg0;
 FyProject_neg_j_1 = FyProject_neg0;    
       Fmax_neg_j_1 = Fmax_neg0;
FmaxProject_neg_j_1 = FmaxProject_neg0;

  beta_S_j_1 = 0;
  beta_C_j_1 = 0;
  beta_K_j_1 = 0;
  beta_F_j_1 = 0;
  
  Energy_Excrsn= 0;
   Energy_Rev= 0;
Energy_total = 0;

    Yield_Flag = 0;
  Fail_FlagPos = 0;
  Fail_FlagNeg = 0;  
Excursion_Flag = 0;
 Reversal_Flag = 0;
    Mrpos_Flag = 0;
    Mrneg_Flag = 0;
   Energy_Flag = 0;  

             Ri_1 = 0;	
             Mi_1 = 0;
             Di_1 = 0;
 Energy_Excrsni_1 = 0;
 
Rreversal=0;
Mreversal=0;
Rintrsct=0;
Mintrsct=0;
RintrsctL=0;
RintrsctR=0;

beta_Sx = 0.0;
Kp_x = Kp_pos_j_1;
Fy_x = Fy_pos_j_1;
Uy_x = Fy_x / Ke;
FyProject_x = Fy_x - Kp_x * Uy_x;


BackbonePos_Original=[Ke; Uy_pos0; Umax_pos0; Kp_pos0; Kpc_pos0; Fy_pos0; FyProject_pos0; Fmax_pos0; FmaxProject_pos0; Uu_pos0; Fres_pos0; Fll_pos0;Ull_pos0;Kul_pos0];
BackboneNeg_Original=[Ke; Uy_neg0; Umax_neg0; Kp_neg0; Kpc_neg0; Fy_neg0; FyProject_neg0; Fmax_neg0; FmaxProject_neg0; Uu_neg0; Fres_neg0; Fll_neg0;Ull_neg0;Kul_neg0];
BackbonePos_j_1=[K_j_1; Uy_pos_j_1; Umax_pos_j_1; Kp_pos_j_1; Kpc_pos_j_1; Fy_pos_j_1; FyProject_pos_j_1; Fmax_pos_j_1; FmaxProject_pos_j_1; Uu_pos0];
BackboneNeg_j_1=[K_j_1; Uy_neg_j_1; Umax_neg_j_1; Kp_neg_j_1; Kpc_neg_j_1; Fy_neg_j_1; FyProject_neg_j_1; Fmax_neg_j_1; FmaxProject_neg_j_1; Uu_neg0];            
Beta_j_1=[beta_S_j_1; beta_C_j_1; beta_K_j_1; beta_F_j_1];
Constants=[Ref_Energy_S; Ref_Energy_C; Ref_Energy_K; c_S; c_C; c_K; D_pos; D_neg];
Flags=[Excursion_Flag; Reversal_Flag; Yield_Flag; Fail_FlagPos; Fail_FlagNeg; Mrpos_Flag; Mrneg_Flag; Energy_Flag];
InitialValues=[Rreversal; Mreversal; Rintrsct; Mintrsct; RintrsctL; RintrsctR; beta_Sx; Kp_x; Fy_x; Uy_x; FyProject_x];
Smooth=Backbone.Smooth;
Smoothparamters=[n; Roffset; LAMBDA_F; c_F; Ref_Energy_F; Smooth];

end