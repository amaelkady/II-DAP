function [fi,ki, Energy,Backbone_j_1,Flag,u0, du_i_1] = IMKPeakOriented_MainFun(Backbone,Energy,Backbone_j_1,Flag,u0,ui, ui_1, du_i_1,fi_1)
%#####################################################################################################################################################################
%IMK PINCHING HYSTERESIS MODEL (Ibarra et al., 2005; Lignos and Krawinkler, 2011)
%#####################################################################################################################################################################
du = ui-ui_1; % Incremental deformation at current step

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
if Flag.Failure ~= 1
    
    %% CHECK FOR UNLOADING
    if     fi_1>0 && du<=0 && du*du_i_1<=0
        Flag.Unload=1;
        Flag.Reversal=1;
        Flag.Reload=0;
        K_check=(Backbone_j_1.FLastPeak_pos-fi_1)/(Backbone_j_1.ULastPeak_pos-ui_1); % In cases where there are small unloading/reloading cycles within an unloading branch
        if K_check >=1.01*Backbone_j_1.Kul || K_check <=0.99*Backbone_j_1.Kul
        Backbone_j_1.FLastPeak_pos=fi_1;
        Backbone_j_1.ULastPeak_pos=ui_1;
        end
    elseif fi_1<0 && du>0 && du*du_i_1<=0
        Flag.Unload=1;
        Flag.Reversal=1;
        Flag.Reload=0;
        K_check=(Backbone_j_1.FLastPeak_neg-fi_1)/(Backbone_j_1.ULastPeak_neg-ui_1);
        if K_check >=1.01*Backbone_j_1.Kul || K_check <=0.99*Backbone_j_1.Kul
        Backbone_j_1.FLastPeak_neg=fi_1;
        Backbone_j_1.ULastPeak_neg=ui_1;
        end
    else
        Flag.Reversal=0;
    end
        
    
    %% CHECK FOR RELOADING
    if     fi_1>0 && du>0 && du_i_1<0
        Flag.Reload=1;
        Flag.Unload=0;
    elseif fi_1<0 && du<0 && du_i_1>0
        Flag.Reload=1;
        Flag.Unload=0;
    end
    
    %% CHECK FOR NEW EXCURSION
    if     fi_1 < 0 && fi_1 + du*Backbone_j_1.Kul >= 0
        Flag.Excursion = 1;
        Flag.Reload=0;
        Flag.Unload=0;
        u0 = ui_1-(fi_1)/Backbone_j_1.Kul; % Deformation at new excursion
    elseif fi_1 > 0 && fi_1 + du*Backbone_j_1.Kul <= 0
        Flag.Excursion = 1;
        Flag.Reload=0;
        Flag.Unload=0;
        u0 = ui_1-(fi_1)/Backbone_j_1.Kul; % Deformation at new excursion
    else
        Flag.Excursion = 0;
    end
    
    
    %% UPDATE GLOBAL PEAK POINTS
    if fi_1>=0 && ui_1 >= Backbone_j_1.Upeak_pos
        Backbone_j_1.Upeak_pos = ui_1;
        Backbone_j_1.Fpeak_pos = fi_1;
    end
    if fi_1<0 &&  ui_1 <= Backbone_j_1.Upeak_neg
        Backbone_j_1.Upeak_neg = ui_1;
        Backbone_j_1.Fpeak_neg = fi_1;
    end
    
    %% CHECK FOR YIELDING
    if     Backbone_j_1.Upeak_pos > Backbone_j_1.Uy_pos || Backbone_j_1.Upeak_neg < Backbone_j_1.Uy_neg
        Flag.Yield = 1;
    end
    
    %% UPDATE DETERIORATION PARAMETERS AT EACH NEW EXCURSION    
    if Flag.Excursion == 1 
        Ei = max(0, Energy.Energy_Acc-Energy.Energy_Diss);                                            % Energy dissipated in previous excursion
        betaS = (Ei/(Backbone.EtS-Energy.Energy_Acc))^Backbone.c_S;
        betaC = (Ei/(Backbone.EtC-Energy.Energy_Acc))^Backbone.c_C;
        betaA = (Ei/(Backbone.EtA-Energy.Energy_Acc))^Backbone.c_A;
        Energy.Energy_Diss=Energy.Energy_Acc;
    else
        betaS = 0;
        betaC = 0;
        betaA = 0;
    end
    
    if Flag.Reversal== 1
        EpjK = Energy.Energy_Acc - 0.5*(fi_1 / Backbone_j_1.Kul)*fi_1;
        EiK = Energy.Energy_Acc - Energy.Energy_Diss + 0.5*(fi_1 / Backbone_j_1.Kul)*fi_1;
        betaK = (EiK / (Backbone.EtK - EpjK))^Backbone.c_K;
        Backbone_j_1.Kul = Backbone_j_1.Kul * (1 - betaK);
    else
        betaK = 0;
    end

    %% UPDATE BACKBONE AND PEAK POINT
    if Flag.Excursion == 1
        % Positive loading backbone
        if fi_1 < 0 && Flag.Yield==1
            % Basic strength deterioration: Yield point
            Backbone_j_1.Uy_pos  = max(Backbone_j_1.Uy_pos -Backbone_j_1.Fy_pos *betaS* Backbone.D_pos/Backbone.Ke,Backbone_j_1.Fres_pos/Backbone.Ke);
            Backbone_j_1.Fy_pos  = max(Backbone_j_1.Fy_pos *(1-betaS* Backbone.D_pos),Backbone_j_1.Fres_pos);
            % Basic strength deterioration: Post-yield Stiffness
            if (Backbone_j_1.Fy_pos ~= Backbone_j_1.Fres_pos)
                Backbone_j_1.Kp_pos = Backbone_j_1.Kp_pos *(1-betaS* Backbone.D_pos);
            else
                Backbone_j_1.Kp_pos  = 0;
            end
            % Basic strength deterioration: Capping Point
            sPCsp = (Backbone_j_1.Fy_pos -Backbone_j_1.Uy_pos *Backbone_j_1.Kp_pos -Backbone_j_1.Fmax_pos +Backbone_j_1.Kpc_pos*Backbone_j_1.Umax_pos )/(Backbone_j_1.Kpc_pos-Backbone_j_1.Kp_pos );
            Backbone_j_1.Fmax_pos  = Backbone_j_1.Fmax_pos +(sPCsp-Backbone_j_1.Umax_pos )*Backbone_j_1.Kpc_pos;
            Backbone_j_1.Umax_pos  = sPCsp;
            % Post-capping strength deterioration: Capping point
            sPCpcp = max(Backbone_j_1.Umax_pos +betaC* Backbone.D_pos*(Backbone_j_1.Fmax_pos -Backbone_j_1.Kpc_pos*Backbone_j_1.Umax_pos )/(Backbone_j_1.Kpc_pos-Backbone_j_1.Kp_pos ),Backbone_j_1.Uy_pos );
            Backbone_j_1.Fmax_pos  = Backbone_j_1.Fmax_pos +(sPCpcp-Backbone_j_1.Umax_pos )*Backbone_j_1.Kp_pos ;
            Backbone_j_1.Umax_pos  = sPCpcp;
            Backbone_j_1.Upeak_pos  = (1+betaA* Backbone.D_pos)*Backbone_j_1.Upeak_pos ;            % Accelerated reloading stiffness deterioration: Target peak deformation point
            if (Backbone_j_1.Upeak_pos  <= Backbone_j_1.Uy_pos )                                    % Target peak deformation in reloading branch of the updated backbone
                Backbone_j_1.Fpeak_pos  = Backbone.Ke*Backbone_j_1.Upeak_pos  ;
            elseif Backbone_j_1.Upeak_pos  <= Backbone_j_1.Umax_pos                                 % Target peak deformation in post-yield branch of the updated backbone
                Backbone_j_1.Fpeak_pos  = Backbone_j_1.Kp_pos *(Backbone_j_1.Upeak_pos -Backbone_j_1.Uy_pos )+Backbone_j_1.Fy_pos ;
            else                                                                                    % Target peak deformation in post-capping branch of the updated backbone
                Backbone_j_1.Fpeak_pos  = max(Backbone_j_1.Kpc_pos*(Backbone_j_1.Upeak_pos -Backbone_j_1.Umax_pos )+Backbone_j_1.Fmax_pos ,Backbone_j_1.Fres_pos);
            end
            
        elseif fi_1 >= 0 && Flag.Yield==1
            
            % Negative loading backbone
            % Basic strength deterioration: Yield point
            Backbone_j_1.Uy_neg  = min(Backbone_j_1.Uy_neg -Backbone_j_1.Fy_neg *betaS* Backbone.D_neg/Backbone.Ke,Backbone_j_1.Fres_neg/Backbone.Ke);
            Backbone_j_1.Fy_neg  = min(Backbone_j_1.Fy_neg *(1-betaS* Backbone.D_neg),Backbone_j_1.Fres_neg);
            % Basic strength deterioration: Post-yield stiffness
            if (Backbone_j_1.Fy_neg  ~= Backbone_j_1.Fres_neg)
                Backbone_j_1.Kp_neg  = Backbone_j_1.Kp_neg *(1-betaS* Backbone.D_neg);
            else
                Backbone_j_1.Kp_neg  = 0;
            end
            % Basic strength deterioration: Capping point
            sPCsn = (Backbone_j_1.Fy_neg -Backbone_j_1.Uy_neg *Backbone_j_1.Kp_neg -Backbone_j_1.Fmax_neg +Backbone_j_1.Kpc_neg*Backbone_j_1.Umax_neg )/(Backbone_j_1.Kpc_neg-Backbone_j_1.Kp_neg );
            Backbone_j_1.Fmax_neg  = Backbone_j_1.Fmax_neg +(sPCsn-Backbone_j_1.Umax_neg )*Backbone_j_1.Kpc_neg;
            Backbone_j_1.Umax_neg  = sPCsn;
            % Post-capping strength deterioration: Capping point
            sPCpcn = min(Backbone_j_1.Umax_neg +betaC* Backbone.D_neg*(Backbone_j_1.Fmax_neg -Backbone_j_1.Kpc_neg*Backbone_j_1.Umax_neg )/(Backbone_j_1.Kpc_neg-Backbone_j_1.Kp_neg ),Backbone_j_1.Uy_neg );
            Backbone_j_1.Fmax_neg  = Backbone_j_1.Fmax_neg +(sPCpcn-Backbone_j_1.Umax_neg )*Backbone_j_1.Kp_neg ;
            Backbone_j_1.Umax_neg  = sPCpcn;
            % Accelerated reloading stiffness deterioration: Target peak deformation point
            Backbone_j_1.Upeak_neg  = (1+betaA* Backbone.D_neg)*Backbone_j_1.Upeak_neg ;
            if (Backbone_j_1.Upeak_neg  >= Backbone_j_1.Uy_neg )                    % Target peak deformation in reloading branch of the updated backbone
                Backbone_j_1.Fpeak_neg  = Backbone.Ke*Backbone_j_1.Upeak_neg ;
            elseif (Backbone_j_1.Upeak_neg  >= Backbone_j_1.Umax_neg )              % Target peak deformation in post-yield branch of the updated backbone
                Backbone_j_1.Fpeak_neg  = Backbone_j_1.Kp_neg *(Backbone_j_1.Upeak_neg -Backbone_j_1.Uy_neg )+Backbone_j_1.Fy_neg ;
            else                                                                    % Target peak deformation in post-capping branch of the updated backbone
                Backbone_j_1.Fpeak_neg  = min(Backbone_j_1.Kpc_neg*(Backbone_j_1.Upeak_neg -Backbone_j_1.Umax_neg )+Backbone_j_1.Fmax_neg ,Backbone_j_1.Fres_neg);
            end

        end
    end
    
    %% Update Deformation at Residual Points
    Backbone_j_1.Ures_pos  = ( Backbone_j_1.Fres_pos - Backbone_j_1.Fmax_pos + Backbone_j_1.Kpc_pos * Backbone_j_1.Umax_pos) / Backbone_j_1.Kpc_pos;
    Backbone_j_1.Ures_neg  = ( Backbone_j_1.Fres_neg - Backbone_j_1.Fmax_neg + Backbone_j_1.Kpc_neg * Backbone_j_1.Umax_neg) / Backbone_j_1.Kpc_neg;

    %% CHECK TARGET POINT: LAST CYCLE PEAK or GLOBAL PEAK
    if Flag.Excursion==1
        if du>=0 
            Krel_LastPeak    = Backbone_j_1.FLastPeak_pos / (Backbone_j_1.ULastPeak_pos - u0);
            Krel_GlobalPeak  = Backbone_j_1.Fpeak_pos     / (Backbone_j_1.Upeak_pos     - u0);
        else
            Krel_LastPeak    = Backbone_j_1.FLastPeak_neg / (Backbone_j_1.ULastPeak_neg - u0);
            Krel_GlobalPeak  = Backbone_j_1.Fpeak_neg     / (Backbone_j_1.Upeak_neg     - u0);
        end
        if du>=0 && Backbone_j_1.FLastPeak_pos >= Backbone_j_1.Fpeak_pos
            Flag.TargetPeak=0;
        elseif du<=0 && Backbone_j_1.FLastPeak_neg <= Backbone_j_1.Fpeak_neg
            Flag.TargetPeak=0;            
        elseif abs(Krel_LastPeak) <= abs(Krel_GlobalPeak)
            Flag.TargetPeak=0;                        
        else
            Flag.TargetPeak=1;
        end
    end
    
    
    %% COMPUTE FORCE INCREMENT
    % Positive Loading
    if fi_1 + du*Backbone_j_1.Kul >= 0
        % CASE 0: At THE ELASTIC SLOPE
        if ui>=0 && Backbone_j_1.Upeak_pos <= Backbone_j_1.Uy_pos && Flag.Yield==0
            if ui >= Backbone_j_1.Uy_pos
                df = Backbone.Ke*(Backbone_j_1.Uy_pos -ui_1) + Backbone_j_1.Kp_pos*(ui -Backbone_j_1.Uy_pos);
            else
                df = du * Backbone.Ke;
            end
            Case=0;

        % CASE 1: EACH NEW EXCURSION
        elseif Flag.Excursion==1
            if Flag.TargetPeak==0
                Backbone_j_1.Krel  =  Backbone_j_1.Fpeak_pos / (Backbone_j_1.Upeak_pos - u0);
            else
                Backbone_j_1.Krel  = Backbone_j_1.FLastPeak_pos / (Backbone_j_1.ULastPeak_pos - u0);
            end
            df = Backbone_j_1.Kul*(u0-ui_1) + Backbone_j_1.Krel*(ui-u0);
            Case=1;

        % CASE 2: WHEN RELOADING
        elseif Flag.Reload==1 && ui <= Backbone_j_1.ULastPeak_pos
            df = du * Backbone_j_1.Kul;
            Case=2;

        % CASE 2: WHEN UNLOADING
        elseif  Flag.Unload==1
            df = du * Backbone_j_1.Kul;
            Case=2;
            
        % CASE 3: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
        elseif Flag.Reload==1 && ui >= Backbone_j_1.ULastPeak_pos && ui <= Backbone_j_1.Upeak_pos && Backbone_j_1.FLastPeak_pos <= Backbone_j_1.Fpeak_pos
            Backbone_j_1.Krel  = (Backbone_j_1.Fpeak_pos-Backbone_j_1.FLastPeak_pos)/(Backbone_j_1.Upeak_pos-Backbone_j_1.ULastPeak_pos);
            if ui_1 <= Backbone_j_1.ULastPeak_pos
                df = Backbone_j_1.Kul*(Backbone_j_1.ULastPeak_pos -ui_1) + Backbone_j_1.Krel*(ui -Backbone_j_1.ULastPeak_pos);
            else
                df = du * Backbone_j_1.Krel;
            end
            Case=3;
    
        % CASE 4: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
        elseif du >= 0 && ((Flag.TargetPeak==0 && ui <= Backbone_j_1.Upeak_pos) || (Flag.TargetPeak==1 && ui <= Backbone_j_1.ULastPeak_pos))
            if Flag.TargetPeak==0
                Backbone_j_1.Krel  = (Backbone_j_1.Fpeak_pos-fi_1)/(Backbone_j_1.Upeak_pos-ui_1);
            else
                Backbone_j_1.Krel  = (Backbone_j_1.FLastPeak_pos-fi_1)/(Backbone_j_1.ULastPeak_pos-ui_1);                
            end
            df = du * Backbone_j_1.Krel;
            Case=4;
        
        % CASE 5: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
        elseif du >= 0 && Flag.TargetPeak==1 && ui >= Backbone_j_1.ULastPeak_pos && ui <= Backbone_j_1.Upeak_pos && Backbone_j_1.FLastPeak_pos <= Backbone_j_1.Fpeak_pos
            
            Backbone_j_1.Krel  = (Backbone_j_1.Fpeak_pos-Backbone_j_1.FLastPeak_pos)/(Backbone_j_1.Upeak_pos-Backbone_j_1.ULastPeak_pos);
            if ui_1 <= Backbone_j_1.ULastPeak_pos
                df = (Backbone_j_1.FLastPeak_pos - fi_1) + Backbone_j_1.Krel*(ui -Backbone_j_1.ULastPeak_pos);
            else
                df = du * Backbone_j_1.Krel;
            end
                
            Case=5;
            
        % CASE 6: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
        elseif du >= 0 && ui <= Backbone_j_1.Umax_pos
            df = du * Backbone_j_1.Kp_pos;
            Case=6;
        
        % CASE 7: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
        elseif du > 0 && ui >= Backbone_j_1.Umax_pos && ui <= Backbone_j_1.Ures_pos
            if ui_1<=Backbone_j_1.Umax_pos && ui>=Backbone_j_1.Umax_pos
                df = Backbone_j_1.Kp_pos*(Backbone_j_1.Umax_pos -ui_1) + Backbone_j_1.Kpc_pos*(ui -Backbone_j_1.Umax_pos);
            else
                df=du * Backbone_j_1.Kpc_pos;
            end
            Case=7;
        
        % CASE 8: WHEN LOADING AND BEYOND THE RESIDUAL POINT
        elseif du > 0 && ui >= Backbone_j_1.Ures_pos
            df = 0.0;
            if Backbone_j_1.Fres_pos == 0
                Flag.Failure = 1;
            end
            Case=8;
        end
        

    end
    
    % Negative Loading
    if fi_1 + du*Backbone_j_1.Kul <= 0
        % CASE 0: At THE ELASTIC SLOPE
        if ui<=0 && Backbone_j_1.Upeak_neg >= Backbone_j_1.Uy_neg && Flag.Yield==0
            if ui <= Backbone_j_1.Uy_neg
                df = Backbone.Ke*(Backbone_j_1.Uy_neg -ui_1) + Backbone_j_1.Kp_neg*(ui -Backbone_j_1.Uy_neg);
            else
                df = du * Backbone.Ke;
            end
            Case=0;

        % CASE 1: EACH NEW EXCURSION
        elseif Flag.Excursion==1
            if Flag.TargetPeak==0
                Backbone_j_1.Krel  = Backbone_j_1.Fpeak_neg / (Backbone_j_1.Upeak_neg - u0);
            else
                Backbone_j_1.Krel  = Backbone_j_1.FLastPeak_neg / (Backbone_j_1.ULastPeak_neg - u0);
            end
            df = Backbone_j_1.Kul*(u0-ui_1) + Backbone_j_1.Krel*(ui-u0);
            Case=1;

        % CASE 2: WHEN RELOADING
        elseif Flag.Reload==1 && ui >= Backbone_j_1.ULastPeak_neg
            df = du * Backbone_j_1.Kul;
            Case=2;

        % CASE 2: WHEN UNLOADING
        elseif Flag.Unload==1
            df = du * Backbone_j_1.Kul;
            Case=2;

        % CASE 3: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
        elseif Flag.Reload==1 && ui <= Backbone_j_1.ULastPeak_neg && ui >= Backbone_j_1.Upeak_neg && Backbone_j_1.FLastPeak_neg >= Backbone_j_1.Fpeak_neg
            Backbone_j_1.Krel  = (Backbone_j_1.Fpeak_neg-Backbone_j_1.FLastPeak_neg)/(Backbone_j_1.Upeak_neg-Backbone_j_1.ULastPeak_neg);
            if ui_1 >= Backbone_j_1.ULastPeak_neg
                df = Backbone_j_1.Kul*(Backbone_j_1.ULastPeak_neg -ui_1) + Backbone_j_1.Krel*(ui -Backbone_j_1.ULastPeak_neg);
            else
                df = du * Backbone_j_1.Krel;
            end
            Case=3;
    
        % CASE 4: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
        elseif du <= 0 && ((Flag.TargetPeak==0 && ui >= Backbone_j_1.Upeak_neg) || (Flag.TargetPeak==1 && ui >= Backbone_j_1.ULastPeak_neg))
            df = du * Backbone_j_1.Krel;
            Case=4;
        
        % CASE 5: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
        elseif du <= 0 && Flag.TargetPeak==1 && ui <= Backbone_j_1.ULastPeak_neg && ui >= Backbone_j_1.Upeak_neg && Backbone_j_1.FLastPeak_neg >= Backbone_j_1.Fpeak_neg
            Backbone_j_1.Krel  = (Backbone_j_1.Fpeak_neg-Backbone_j_1.FLastPeak_neg)/(Backbone_j_1.Upeak_neg-Backbone_j_1.ULastPeak_neg);
            if ui_1 >= Backbone_j_1.ULastPeak_neg
                df = (Backbone_j_1.FLastPeak_neg - fi_1) + Backbone_j_1.Krel*(ui -Backbone_j_1.ULastPeak_neg);
            else
                df = du * Backbone_j_1.Krel;
            end
            Case=5;
            
        % CASE 6: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
        elseif du <= 0 && ui >= Backbone_j_1.Umax_neg
            df = du * Backbone_j_1.Kp_neg;
            Case=6;
        
        % CASE 7: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
        elseif du < 0 && ui <= Backbone_j_1.Umax_neg && ui >= Backbone_j_1.Ures_neg
            if ui_1>=Backbone_j_1.Umax_neg && ui<=Backbone_j_1.Umax_neg
                df = Backbone_j_1.Kp_neg*(Backbone_j_1.Umax_neg -ui_1) + Backbone_j_1.Kpc_neg*(ui -Backbone_j_1.Umax_neg);
            else
                df=du * Backbone_j_1.Kpc_neg;
            end
            Case=7;
        
        % CASE 8: WHEN LOADING AND BEYOND THE RESIDUAL POINT
        elseif du < 0 &&  ui <= Backbone_j_1.Ures_neg
            df = 0.0;
            if Backbone_j_1.Fres_neg == 0
                Flag.Failure = 1;
            end
            Case=8;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% Force
    fi = fi_1+df;
    
    %% Failure
    % Failure criteria (Tolerance = 1%)
    FailS = (betaS < -0.01 || betaS > 1.01);
    FailC = (betaC < -0.01 || betaC > 1.01);
    FailA = (betaA < -0.01 || betaA > 1.01);
    FailK = (betaK < -0.01 || betaK > 1.01);
    if (FailS || FailC || FailA || FailK)
        fi = 0;
        Flag.Failure = 1;
    end
    if ui>=0.0 && ui >= Backbone.Uu_pos
        fi = 0;
        Flag.Failure=1;
    elseif ui<0.0 && ui <= -Backbone.Uu_neg
        fi = 0;
        Flag.Failure=1;
    end
    if Backbone_j_1.Fpeak_pos==0 || Backbone_j_1.Fpeak_neg==0
        fi=0;
        Flag.Failure=1;
    end
    
    dEi = 0.5*(fi+fi_1)*du; % Internal energy increment

else
    
    fi = 0; % Force at failure
    dEi = 0;
    
end


%% Energies
Energy.Energy_Acc  = Energy.Energy_Acc+dEi;   % Total internal energy accumulated until current increment

%% Update Variables
du_i_1=du;

if fi==fi_1
    ki=10^-6;
elseif du==0
    ki=Backbone.Ke;
else
    ki=(fi-fi_1)/(du);
end
