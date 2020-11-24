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

function [Ri, Mi, Di, Di_1, BackbonePos_j_1, BackboneNeg_j_1, Beta_j_1, Energy_Excrsni_1, Energy_Excrsn, Energy_Rev, Energy_total, Flags, InitialValues, Smoothparamters, TangentK]=IMKSmooth_MainFun (Ri, Ri_1, Mi_1, Di_1, BackbonePos_Original, BackboneNeg_Original, BackbonePos_j_1, BackboneNeg_j_1 , Beta_j_1, Energy_Excrsni_1, Energy_Excrsn, Energy_Rev, Energy_total, Flags, Constants, InitialValues, Smoothparamters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  UNPACK  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

              Ke = BackbonePos_Original(1);
   Kpc_pos0 = BackbonePos_Original(5);
        Fy_pos0 = BackbonePos_Original(6);
         Fres_pos0 = BackbonePos_Original(11);
    if Fres_pos0<0.0; Fres_pos0=0; end
         
   Kpc_neg0 = BackboneNeg_Original(5);
        Fy_neg0 = BackboneNeg_Original(6);
         Fres_neg0 = BackboneNeg_Original(11);
    if Fres_neg0<0.0; Fres_neg0=0; end

              K_j_1 = BackbonePos_j_1(1);
    Uy_pos_j_1 = BackbonePos_j_1(2);
  Umax_pos_j_1 = BackbonePos_j_1(3);
    Kp_pos_j_1 = BackbonePos_j_1(4);
   Kpc_pos_j_1 = BackbonePos_j_1(5);
        Fy_pos_j_1 = BackbonePos_j_1(6);
 FyProject_pos_j_1 = BackbonePos_j_1(7);
       Fmax_pos_j_1 = BackbonePos_j_1(8);
FmaxProject_pos_j_1 = BackbonePos_j_1(9);
        Uu_pos = BackbonePos_j_1(10);

    Uy_neg_j_1 = BackboneNeg_j_1(2);
  Umax_neg_j_1 = BackboneNeg_j_1(3);
    Kp_neg_j_1 = BackboneNeg_j_1(4);
   Kpc_neg_j_1 = BackboneNeg_j_1(5);
        Fy_neg_j_1 = BackboneNeg_j_1(6);
 FyProject_neg_j_1 = BackboneNeg_j_1(7);
       Fmax_neg_j_1 = BackboneNeg_j_1(8);
FmaxProject_neg_j_1 = BackboneNeg_j_1(9);
        Uu_neg = BackboneNeg_j_1(10);

  beta_S_j_1 = Beta_j_1(1);
  beta_C_j_1 = Beta_j_1(2);
  beta_K_j_1 = Beta_j_1(3);
  beta_F_j_1 = Beta_j_1(4);
  
Ref_Energy_S = Constants(1);
Ref_Energy_C = Constants(2);
Ref_Energy_K = Constants(3);
       c_S = Constants(4);
       c_C = Constants(5);
       c_K = Constants(6);
     D_pos = Constants(7);
     D_neg = Constants(8);

Excursion_Flag = Flags(1);
 Reversal_Flag = Flags(2);
    Yield_Flag = Flags(3);
  Fail_FlagPos = Flags(4);
  Fail_FlagNeg = Flags(5);
    Mrpos_Flag = Flags(6);  
    Mrneg_Flag = Flags(7);  
    Energy_Flag = Flags(8);  

Rreversal = InitialValues(1);
Mreversal = InitialValues(2);
 Rintrsct = InitialValues(3);
 Mintrsct = InitialValues(4);
RintrsctL = InitialValues(5);
RintrsctR = InitialValues(6);
beta_Sx = InitialValues(7);
 Kp_x = InitialValues(8);
 Fy_x = InitialValues(9);
Uy_x = InitialValues(10);
FyProject_x = InitialValues(11);

 n = Smoothparamters(1);
 Roffset = Smoothparamters(2);
LAMBDA_F = Smoothparamters(3);
c_F = Smoothparamters(4);
Ref_Energy_F= Smoothparamters(5);
Smooth= Smoothparamters(6);

Mi_boundary_Pos=0;
Mi_boundary_Neg=0;
QuarterFlag=0;
Reversal_Flag=0;


Fmax_pos_j= Fmax_pos_j_1;
Fy_pos_j=Fy_pos_j_1;
FyProject_pos_j=FyProject_pos_j_1;
FmaxProject_pos_j=FmaxProject_pos_j_1;
Kp_pos_j=Kp_pos_j_1;
Kpc_pos_j=Kpc_pos_j_1;
Uy_pos_j=Uy_pos_j_1;		
Umax_pos_j=Umax_pos_j_1;
         
Fy_neg_j=Fy_neg_j_1;
FyProject_neg_j=FyProject_neg_j_1;
FmaxProject_neg_j=FmaxProject_neg_j_1;
Fmax_neg_j=Fmax_neg_j_1;
Kp_neg_j=Kp_neg_j_1;
Kpc_neg_j=Kpc_neg_j_1;
Uy_neg_j=Uy_neg_j_1;		
Umax_neg_j=Umax_neg_j_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find the direction of the current increment "Di": 1:if Ri is moving right    and -1: if Ri is moving left
    if Ri>=Ri_1
        Di=1;
	else
		Di=-1;
    end
    
    % Simple Notation for current step parameters
	Mi= Mi_1 + K_j_1*(Ri-Ri_1);    
    
    % Check for Fail Flag
    if Ri>=Uu_pos
        Fail_FlagPos=1;
    end
    if Ri<=-Uu_neg
        Fail_FlagNeg=1;
    end
    
    % Get Information before first Yield
    if Mi > Fy_pos0  && Yield_Flag==0
       Yield_Flag=1;
    end

	% Check if previous point was a reversal point
	if Di_1/Di<0
		Reversal_Flag=1;
        Rreversal=Ri_1;
        Mreversal=Mi_1;
    end

    if Reversal_Flag==1
		Rintrsct_K=Rreversal-Mreversal/K_j_1;
        Energy_Rev=max(0,Energy_total-Energy_Excrsni_1+0.5*Mreversal*(Rintrsct_K-Rreversal)); % total energy dissipated right before this drift reversal (minus the elastic energy)
        % Update Loading/Unloading Stiffness
		beta_K_j   =(Energy_Rev/(2*Ref_Energy_K-Energy_total+0.5*Mreversal*(Rintrsct_K-Rreversal)))^c_K;
        K_j=K_j_1*(1-beta_K_j);
        if Mrpos_Flag==1 || Mrneg_Flag==1
            K_j=0.5*Ke;
        end
        % Update Paramters for Smooth Transition Option
		beta_F_j   =(Energy_Rev/(2*Ref_Energy_F-Energy_total+0.5*Mreversal*(Rintrsct_K-Rreversal)))^c_F;
		if Di>0
            beta_Sx = (max(0,(Energy_Rev + 0.5*Mreversal* (Rintrsct_K - Rreversal)))/ (Ref_Energy_S - Energy_total))^ c_S;
			Kp_x = Kp_pos_j_1   * (1.0 - beta_Sx * D_pos);
			Fy_x = Fy_pos_j_1 	* (1.0 - beta_Sx * D_pos);
			Uy_x = Fy_x / K_j;
			FyProject_x = Fy_x - Kp_x * Uy_x;
        else
            beta_Sx = (max(0,(Energy_Rev + 0.5*Mreversal* (Rintrsct_K - Rreversal)))/ (Ref_Energy_S - Energy_total))^ c_S;
			Kp_x = Kp_neg_j_1   * (1.0 - beta_Sx * D_neg);
			Fy_x = Fy_neg_j_1 	* (1.0 - beta_Sx * D_neg);
			Uy_x = Fy_x / K_j;
			FyProject_x = Fy_x - Kp_x * Uy_x;
        end
        
    else
        beta_F_j = beta_F_j_1;
        beta_K_j = beta_K_j_1;
		K_j=K_j_1;
 	end
    %%
    % Calculate Backbone paramters at current excursion based on Energy Disspipated in the previous Excursion 
    if Excursion_Flag==1
        beta_S_j =(Energy_Excrsn/(Ref_Energy_S-Energy_total))^c_S;
        beta_C_j =(Energy_Excrsn/(Ref_Energy_C-Energy_total))^c_C;
      
        if Ri-Ri_1 >=0.0
            % Update Fy, Fmax Projection, Kp, and Kpc for Current Step
            Fy_pos_j= Fy_pos_j_1  * (1-beta_S_j*D_pos);
            FmaxProject_pos_j=FmaxProject_pos_j_1 * (1-beta_C_j*D_pos);
            Kp_pos_j=Kp_pos_j_1* (1-beta_S_j*D_pos);
            if  Fres_pos0==0
                Kpc_pos_j=Kpc_pos0 * (FmaxProject_pos_j-Fres_pos0)/FmaxProject_pos_j;               
            else
                Kpc_pos_j=Kpc_pos0 * (Fy_pos_j-Fres_pos0)/(Fy_pos0-Fres_pos0);        
            end
            % Calculate Rotation at Capping Point (Intersection of the two Slopes)        
            Uy_pos_j=Fy_pos_j/K_j;
            FyProject_pos_j=Fy_pos_j-Kp_pos_j*Uy_pos_j;
            Umax_pos_j=abs((FmaxProject_pos_j - FyProject_pos_j) / (Kpc_pos_j+Kp_pos_j));
            Fmax_pos_j=FyProject_pos_j+Umax_pos_j*Kp_pos_j;

        else
            % Update Fy, Fmax Projection, Kp, and Kpc for Current Step
            Fy_neg_j= Fy_neg_j_1 * (1-beta_S_j*D_neg);
            FmaxProject_neg_j=FmaxProject_neg_j_1 * (1-beta_C_j*D_neg);
            Kp_neg_j=Kp_neg_j_1* (1-beta_S_j*D_neg);
            if  Fres_pos0==0
                Kpc_neg_j=Kpc_neg0 * (FmaxProject_neg_j-Fres_neg0)/FmaxProject_neg_j;
            else
                Kpc_neg_j=Kpc_neg0 * (Fy_neg_j-Fres_neg0)/(Fy_neg0-Fres_neg0);
            end
            % Calculate Rotation at Capping Point (Intersection of the two Slopes)        
            Uy_neg_j=Fy_neg_j/K_j;
            FyProject_neg_j=Fy_neg_j-Kp_neg_j*Uy_neg_j;
            Umax_neg_j=abs((FmaxProject_neg_j - FyProject_neg_j) / (Kpc_neg_j+Kp_neg_j));
            Fmax_neg_j=FyProject_neg_j+Umax_neg_j*Kp_neg_j;            
        end
    
    else
        beta_S_j =beta_S_j_1;
        beta_C_j=beta_C_j_1;

		if Di>=0.0
			Fy_pos_j=Fy_pos_j_1;
			FyProject_pos_j=FyProject_pos_j_1;
			FmaxProject_pos_j=FmaxProject_pos_j_1;
			Fmax_pos_j=Fmax_pos_j_1;
			Kp_pos_j=Kp_pos_j_1;
			Kpc_pos_j=Kpc_pos_j_1;
			Uy_pos_j=Uy_pos_j_1;		
			Umax_pos_j=Umax_pos_j_1;
		else
			Fy_neg_j=Fy_neg_j_1;
			FyProject_neg_j=FyProject_neg_j_1;
			FmaxProject_neg_j=FmaxProject_neg_j_1;
			Fmax_neg_j=Fmax_neg_j_1;
			Kp_neg_j=Kp_neg_j_1;
			Kpc_neg_j=Kpc_neg_j_1;
			Uy_neg_j=Uy_neg_j_1;		
			Umax_neg_j=Umax_neg_j_1;
		end
    end

    %% If the residual moment is reached in a given direction, Override the values of Fmax, Umax, Kp and Kpc
    if Di>=0.0
        if Fmax_pos_j <Fres_pos0
           Fmax_pos_j= Fres_pos0;          
           Kpc_pos_j=10^-6;
           Kp_pos_j=10^-6;
		   Umax_pos_j=10^-6;
        end      
	else
        if Fmax_neg_j <Fres_neg0        
           Fmax_neg_j= Fres_neg0;
           Kpc_neg_j=10^-6;
           Kp_neg_j=10^-6;
           Umax_neg_j=10^-6;
        end
    end
    
    % Simple and unified notation for current bacbone parameters
	if Di >=0.0
		Ki=K_j;
		Kpi=Kp_pos_j;
		Kpci=Kpc_pos_j;
		Uyi=Uy_pos_j;
		Umaxi=Umax_pos_j;
		Fmaxi=Fmax_pos_j;
		FyProjecti=FyProject_pos_j;
		FmaxProjecti=FmaxProject_pos_j;
		Upci=Fmaxi/Kpci;
	else
		Ki=K_j;
		Kpi=Kp_neg_j;
		Kpci=Kpc_neg_j;
		Uyi=Uy_neg_j;
		Umaxi=Umax_neg_j;
		Fmaxi=Fmax_neg_j;
		FyProjecti=FyProject_neg_j;
		FmaxProjecti=FmaxProject_neg_j;
		Upci=Fmaxi/Kpci;
    end

	
    % Moment Calculation Based on unloading/reloading stiffeness
    Mi= Mi_1 + Ki*(Ri-Ri_1);
    
    %% Location Flags
    if Ri>=0 && Mi>0
        QuarterFlag=1;
    elseif Ri>=0 && Mi<0
        QuarterFlag=2;
    elseif Ri<=0 && Mi<0
        QuarterFlag=3;
    elseif Ri<=0 && Mi>0
        QuarterFlag=4;
    end
    
    % Get Boundary Moment at Current Step Based on Current BackBone Curve
    if Ri >= 0.0 &&  abs(Ri)<=Umaxi
        Mi_boundary_Pos=  FyProjecti+Kpi*Ri;
        Mi_boundary_Neg= -FyProjecti+Kpi*Ri;
    elseif Ri >= 0.0 && abs(Ri)>Umaxi
        Mi_boundary_Pos=  FmaxProjecti-Kpci*Ri;
        Mi_boundary_Neg= -FyProjecti+Kpi*Ri;     
    elseif Ri <= 0.0 && abs(Ri)<=Umaxi 
        Mi_boundary_Pos=  FyProjecti+Kpi*Ri;
        Mi_boundary_Neg= -FyProjecti+Kpi*Ri;
    elseif Ri <= 0.0 && abs(Ri)>Umaxi
        Mi_boundary_Pos=  FyProjecti+Kpi*Ri;     
        Mi_boundary_Neg= -FmaxProjecti-Kpci*Ri;    
    end

   
    %% Check and Modify Boundary Moment if it Exceeds Mr
    if QuarterFlag==3 && Mi_boundary_Neg >= 0.0
      Mi_boundary_Neg= 0.0;
    end
    if QuarterFlag==3 && abs(Mi_boundary_Neg) <= Fres_neg0 
      Mi_boundary_Neg= -Fres_neg0;
      Mrneg_Flag=1;
    end 
    if QuarterFlag==1 && Mi_boundary_Pos <= 0.0 
      Mi_boundary_Pos= 0.0;
    end
    if QuarterFlag==1 && abs(Mi_boundary_Pos) <= Fres_pos0
      Mi_boundary_Pos= Fres_pos0;
      Mrpos_Flag=1;
    end

    
    
Mi_boundary_Pos1=99999999;
Mi_boundary_Neg1=-99999999;    
if Smooth==1
%%////////////////////////////////////////////////////////////////////
%//////////////////// CODE FOR SMOOTH SPLINE //////////////////////////
%//////////////////////////////////////////////////////////////////////
%/  Get Boundary Moment in the Transition Region at Current Step Based on Spline Curve

if Di < 0
    b1 = Mreversal - Ki*Rreversal;
    b2 = -FyProject_x;
    RintrsctXaxis = -b1 / Ki;
    Rintrsct = -(b2 - b1) / (Kp_x - Ki);

    RintrsctL = Rintrsct - Roffset / (1 - beta_F_j);
    if RintrsctL > Umaxi
        RintrsctL = Umaxi;
    end
    MintrsctL = b2 + Kpi*RintrsctL;

    RintrsctR = Rintrsct + n*Roffset / (1 - beta_F_j);
    if RintrsctR > Rreversal
        RintrsctR = Rreversal;
    end
    MintrsctR = b1 + Ki*RintrsctR;

    if Ri>RintrsctL && Ri<RintrsctR
        % R and M vectors for PCHIP interpolation
        Rx = [RintrsctL-0.01           RintrsctL   RintrsctR           RintrsctR+0.01];
        Mx = [b2+Kpi*(RintrsctL-0.01)  MintrsctL   MintrsctR    b1+Ki*(RintrsctR+0.01)];
        Mi_boundary_Pos1 =interp1(Rx,Mx,Ri,'pchip');
        Mi_boundary_Neg1 =interp1(Rx,Mx,Ri,'pchip');
    end
end

if Di >= 0
    b1 = Mreversal - Ki*Rreversal;
    b2 = FyProject_x;
    RintrsctXaxis = -b1 / Ki;
    Rintrsct = -(b2 - b1) / (Kp_x - Ki);

    RintrsctL = Rintrsct - n*Roffset / (1 - beta_F_j);
    if RintrsctL <Rreversal
        RintrsctL = Rreversal;
    end
    MintrsctL = b1 + Ki*RintrsctL;

    RintrsctR = Rintrsct + Roffset / (1 - beta_F_j);
    if RintrsctR >Umaxi
        RintrsctR = Umaxi;
    end
    MintrsctR = b2 + Kpi*RintrsctR;

    if Ri>RintrsctL && Ri<RintrsctR
        % R and M vectors for PCHIP interpolation
        Rx= [RintrsctL-0.01           RintrsctL    RintrsctR            RintrsctR+0.01];
        Mx= [b1+Ki*(RintrsctL-0.01)   MintrsctL    MintrsctR    b2+Kpi*(RintrsctR+0.01)];
        Mi_boundary_Pos1 = interp1(Rx,Mx,Ri,'pchip');
        Mi_boundary_Neg1 = interp1(Rx,Mx,Ri,'pchip');
    end
end

end
%%///////////////////////////////////////////////////////////////////////
%%///////////////////////////////////////////////////////////////////////
if Mi_boundary_Pos>Mi_boundary_Pos1; Mi_boundary_Pos=Mi_boundary_Pos1; end
if Mi_boundary_Neg<Mi_boundary_Neg1; Mi_boundary_Neg=Mi_boundary_Neg1; end

%         h40=scatter(RintrsctL,MintrsctL,'og');
%         h41=scatter(RintrsctR,MintrsctR,'^g');
%         delete (h40)
%         delete (h41)
	% If Failure took place in a given direction (Fail_Flag_dir=1), Set the Boundary Moment in the opposite direction to Mr
    if Ri<0 && Di>0 && Fail_FlagNeg==1
      Mi_boundary_Pos= Fres_pos0;
    end
    if Ri>0 && Di<0 && Fail_FlagPos==1
      Mi_boundary_Neg= -Fres_neg0;
    end

	% If the residual moment is reached in a given direction (Mrdir_Flag=1), Set the Boundary Moment in both directions to Mr
    if Mrpos_Flag==1 || Mrneg_Flag == 1
      Mi_boundary_Pos=  Fres_pos0;
      Mi_boundary_Neg= -Fres_neg0;
    end
    
	%%%%%%% Current Step Moment Calculation %%%%%%%
	% If current moment based on unloading/reloading Ki is larger than the boundary moment, set it equal to the boundary moment
	if  Mi > Mi_boundary_Pos
        Mi = Mi_boundary_Pos;
    end
    if  Mi < Mi_boundary_Neg
        Mi = Mi_boundary_Neg; 
    end 

    % if fail flag is reached in a given loading direction, set current moment equal to zero 
    if Ri>0 && Fail_FlagPos==1
       Mi=0.0;
    end
    if Ri<0 && Fail_FlagNeg==1
       Mi=0.0;
    end    

    if Energy_Flag==1
       Mi=0.0;
    end 
    
    %% Energy Calculation
    Energy_total=Energy_total+(Mi+Mi_1)*0.5*(Ri-Ri_1); % total energy dissipated till current incremental step
    
    if Mi/Mi_1<0
        Energy_Excrsn=max(0,Energy_total-Energy_Excrsni_1);  % total energy dissipated in current excursion
        Energy_Excrsni_1=Energy_total;      		  % total energy dissipated in previous excursion
        Excursion_Flag=1;
    else
        Excursion_Flag=0;
    end

	% Check if the Component inheret Reference Energy is Consumed
    if Excursion_Flag==1
        if Energy_total >= min ([Ref_Energy_S; Ref_Energy_C])
           Energy_Flag=1; 
        end
        if beta_S_j>1 || beta_C_j>1
           Energy_Flag=1;         
        end
    elseif Reversal_Flag==1
        if Energy_total >= Ref_Energy_K
           Energy_Flag=1; 
        end      
        if beta_K_j>1
           Energy_Flag=1;         
        end
        if Smooth==1
            if Energy_total >= Ref_Energy_F
               Energy_Flag=1; 
            end      
            if beta_F_j>1
               Energy_Flag=1;         
            end
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RETURN VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Di_1=Di;
if Ri>=Ri_1
BackbonePos_j_1=[K_j; Uy_pos_j; Umax_pos_j; Kp_pos_j; Kpc_pos_j; Fy_pos_j; FyProject_pos_j; Fmax_pos_j; FmaxProject_pos_j; Uu_pos];
else
BackboneNeg_j_1=[K_j; Uy_neg_j; Umax_neg_j; Kp_neg_j; Kpc_neg_j; Fy_neg_j; FyProject_neg_j; Fmax_neg_j; FmaxProject_neg_j; Uu_neg];            
end
Beta_j_1=[beta_S_j; beta_C_j; beta_K_j; beta_F_j];
Flags=[Excursion_Flag; Reversal_Flag; Yield_Flag; Fail_FlagPos; Fail_FlagNeg; Mrpos_Flag; Mrneg_Flag; Energy_Flag];
InitialValues=[Rreversal; Mreversal; Rintrsct; Mintrsct; RintrsctL; RintrsctR; beta_Sx; Kp_x; Fy_x; Uy_x; FyProject_x];
Smoothparamters=[n; Roffset; LAMBDA_F; c_F; Ref_Energy_F; Smooth];
%Others=[Mi_boundary_Pos; Mi_boundary_Neg; QuarterFlag; Reversal_Flag; FyProject_x; Kp_x];
% Tangent Stiffeness Calculation
if Ri == Ri_1
	TangentK = 0.0;
	Mi = Mi_1;
else
	TangentK = (Mi - Mi_1) / (Ri - Ri_1);
end

if Mi == Fres_pos0
	TangentK = 10^-6;
end