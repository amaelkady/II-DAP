function [Teq,U,V,A,F,K,Ag,Flag_PDeltaFail,PDeltaFail_Index] = Fun_ResponseHistory (GM_No,Sa_SF, MassT,TT,CT)
global MainDirectory ProjectPath ProjectName

cd (ProjectPath)
load (ProjectName)
load (ProjectName,'gamma','beta')
cd (MainDirectory)
    
Flag_PDeltaFail=0;
M=MassT;
T=TT;
C=CT;

nRefinePoints=GM.dt(GM_No)/dt_user;

evalc([' Ag1(:,1)=GM.GA',  num2str(GM_No),'*Sa_SF']);
evalc(['Teq1(:,1)=GM.Time',num2str(GM_No)]);

timex=0:1/(size(Ag1,1)-1):1;  timex=timex';
x=0:1/(size(timex,1)*nRefinePoints-1):1;  x=x';
Ag =pchip(timex,Ag1, x);
Teq=pchip(timex,Teq1,x);

Teq(1:end+1,1)=[0.0; Teq(:,1)];
Ag (1:end+1,1)=[0.0;  Ag(:,1)];

PDeltaFail_Index=size(Teq,1)+1;

ui_1=0;
ui=0;
vi=0;
ai=0;
ki=Backbone.Ke;
fi=0;

U(1,1)=0;
V(1,1)=0;
A(1,1)=0;
F(1,1)=0;
K(1,1)=Backbone.Ke;
% 
% figure
% hold on; grid on; box on;
% set(gca,'Ylim',[-2 2]); 
% set(gca,'Xlim',[-0.5 0.5]);
%     UFail_pos0=min(Backbone.Uu_pos0,Backbone.Uy_pos0+Backbone.Up_pos0+Backbone.Upc_pos0);
%     UFail_neg0=min(Backbone.Uu_neg0,Backbone.Uy_neg0+Backbone.Up_neg0+Backbone.Upc_neg0);
% 
%     % Plot backbone
%     if Backbone.Uu_pos0 > Backbone.Uy_pos0+Backbone.Up_pos0+Backbone.Upc_pos0
%         RotationPos=[0 Backbone.Uy_pos0 Backbone.Umax_pos0 Backbone.Umax_pos0+(Backbone.Fmax_pos0-Backbone.Fres_pos0)/Backbone.Kpc_pos0 Backbone.Uu_pos0 Backbone.Uu_pos0];
%         MomentPos=[0 Backbone.Fy_pos0 Backbone.Fmax_pos0 Backbone.Fres_pos0 Backbone.Fres_pos0 0];
%     else
%         RotationPos=[0 Backbone.Uy_pos0 Backbone.Umax_pos0 Backbone.Uu_pos0 Backbone.Uu_pos0];
%         MomentPos  =[0 Backbone.Fy_pos0 Backbone.Fmax_pos0 Backbone.Fmax_pos0-Backbone.Kpc_pos0*(Backbone.Uu_pos0-Backbone.Umax_pos0)      0];    
%     end
% 
%     if Backbone.Uu_neg0 > Backbone.Uy_neg0+Backbone.Up_neg0+Backbone.Upc_neg0
%         RotationNeg=[0 Backbone.Uy_neg0 Backbone.Umax_neg0 Backbone.Umax_neg0+(Backbone.Fmax_neg0-Backbone.Fres_neg0)/Backbone.Kpc_neg0 Backbone.Uu_neg0 Backbone.Uu_neg0];
%         MomentNeg=[0 Backbone.Fy_neg0 Backbone.Fmax_neg0 Backbone.Fres_neg0 Backbone.Fres_neg0 0];
%     else
%         RotationNeg=[0 Backbone.Uy_neg0 Backbone.Umax_neg0 Backbone.Uu_neg0 Backbone.Uu_neg0];
%         MomentNeg  =[0 Backbone.Fy_neg0     Backbone.Fmax_neg0      Backbone.Fmax_neg0-Backbone.Kpc_neg0*(Backbone.Uu_neg0-Backbone.Umax_neg0)      0];    
%     end
%     
%     % Plot backbone
%     plot(RotationPos,MomentPos,'-k','linewidth',2);
%     plot(-RotationNeg,-MomentNeg,'-k','linewidth',2);
%     
for i=1:size(Ag,1)-1

    %% Time step and external forces (ground accelerations) at start and end of current step
    agi =Ag(i  ,1);
    agii=Ag(i+1,1);
    dt = Teq(i+1,1)-Teq(i,1);

    %% Run NewMark
    if Integrator_Option==1
        [uii,vii,aii, ui_1] = Integrator_CentralDifferenceSingleIncrNL(C,M,dt,ui_1,ui,agii,fi);
    elseif Integrator_Option==2   
        [uii,vii,aii] = Integrator_NewmarkSingleIncrNL(ki,C,M,gamma,beta,dt,ui,vi,fi,agi,agii);
    elseif Integrator_Option==3 
        [uii,vii,aii] = Integrator_ReccuranceFormulasSingleIncr(ki,T,M,zeta,dt,ui,vi,agi,agii);
    end
    %% Get force at end of current step
    if SystemModel_Option==1
        [kii,fii]=System_Linear (Backbone.Ke, uii,Backbone.du);
    elseif SystemModel_Option==2
        [kii,fii]=System_BiLinear (Backbone.Ke, Backbone.dy, Backbone.Kp, fi, ui, uii);
    elseif SystemModel_Option==3
        % initializing function
        if i==1
            [BackbonePos_Original, BackboneNeg_Original, BackbonePos_j_1, BackboneNeg_j_1, Beta_j_1, Constants, Flags, InitialValues, Smoothparamters, ui, fi, Di, Energy_Excrsni_1,Energy_Excrsn, Energy_Rev, Energy_total]=IMKSmooth_IntializeFun (Backbone);
        end
        % main function
        [uii, fii, Dii, Di , BackbonePos_j_1, BackboneNeg_j_1, Beta_j_1,  Energy_Excrsni_1,Energy_Excrsn, Energy_Rev, Energy_total, Flags, InitialValues, Smoothparamters, kii]=IMKSmooth_MainFun (uii, ui, fi, Di, BackbonePos_Original, BackboneNeg_Original, BackbonePos_j_1, BackboneNeg_j_1 , Beta_j_1, Energy_Excrsni_1,Energy_Excrsn, Energy_Rev, Energy_total, Flags, Constants, InitialValues, Smoothparamters);
    elseif SystemModel_Option==4
        % initializing function
        if i==1
            [Energy,Backbone_j_1,Flag,u0, du_i_1] = IMKPinching_IntializeFun(Backbone);
        end
        % main function
        [fii,kii, Energy,Backbone_j_1,Flag,u0, du_i_1] = IMKPinching_MainFun(Backbone,Energy,Backbone_j_1,Flag,u0,uii, ui, du_i_1,fi);
    elseif SystemModel_Option==5
        % initializing function
        if i==1
            [Energy,Backbone_j_1,Flag,u0, du_i_1] = IMKPeakOriented_IntializeFun(Backbone);
        end
        % main function
        [fii,kii, Energy,Backbone_j_1,Flag,u0, du_i_1] = IMKPeakOriented_MainFun(Backbone,Energy,Backbone_j_1,Flag,u0,uii, ui, du_i_1,fi);
        
    elseif SystemModel_Option==6
        % initializing function
        if i==1
            [BackbonePos_Original, BackboneNeg_Original, BackbonePos_j_1, BackboneNeg_j_1, Beta_j_1, Constants, Flags, InitialValues, Smoothparamters, ui, fi, Di, Energy_Excrsni_1,Energy_Excrsn, Energy_Rev, Energy_total]=IMKSelfCentering_IntializeFun (Backbone);
        end
        % main function
        [uii, fii, Dii, Di , BackbonePos_j_1, BackboneNeg_j_1, Beta_j_1,  Energy_Excrsni_1,Energy_Excrsn, Energy_Rev, Energy_total, Flags, InitialValues, Smoothparamters, kii]=IMKSelfCentering_MainFun (uii, ui, fi, Di, BackbonePos_Original, BackboneNeg_Original, BackbonePos_j_1, BackboneNeg_j_1 , Beta_j_1, Energy_Excrsni_1,Energy_Excrsn, Energy_Rev, Energy_total, Flags, Constants, InitialValues, Smoothparamters);
   end
    
%         h1= plot([0  BackbonePos_j_1(2)  BackbonePos_j_1(3)  BackbonePos_j_1(3)+(BackbonePos_j_1(8)-BackbonePos_Original(11))/BackbonePos_j_1(5) ],[0  BackbonePos_j_1(6)  BackbonePos_j_1(8)  BackbonePos_Original(11) ],'-ok');
%         h2= plot([0 -BackboneNeg_j_1(2) -BackboneNeg_j_1(3) -BackboneNeg_j_1(3)-(BackboneNeg_j_1(8)-BackboneNeg_Original(11))/BackboneNeg_j_1(5) ],[0 -BackboneNeg_j_1(6) -BackboneNeg_j_1(8) -BackboneNeg_Original(11) ],'-ok');
%         plot([Backbone.Uy_pos Backbone.Uy_pos],[-2 2],'--b');
%         plot([Backbone.Uy_neg Backbone.Uy_neg],[-2 2],'--b');
%         h3=plot(U(1:i,1),F(1:i,1),'-r');
%         h4=scatter(U(i,1),F(i,1),'ob');
%         pause(0.01)
%         delete (h1)
%         delete (h2)
%         delete (h3)
%         delete (h4)


    %% Check for Failure due to PDelta
    if PDeltaFail_Index==size(Teq,1)+1
        if SystemModel_Option==1 || SystemModel_Option==2
            if abs(uii) >= delta_PDeltaFail || abs(uii) >= delta_SystemFail
                Flag_PDeltaFail=1;
                PDeltaFail_Index=i;
            end
        elseif SystemModel_Option==3 || SystemModel_Option==4 || SystemModel_Option==5 || SystemModel_Option==6
            if uii>0
                if abs(uii) >= delta_PDeltaFail_pos || abs(uii) >= delta_SystemFail_pos
                    Flag_PDeltaFail=1;
                    delta_PDeltaFail=delta_PDeltaFail_pos;
                    delta_SystemFail=delta_SystemFail_pos;
                    PDeltaFail_Index=i;
                end
            else
                if abs(uii) >= delta_PDeltaFail_neg || abs(uii) >= delta_SystemFail_neg
                    Flag_PDeltaFail=1;
                    delta_PDeltaFail=delta_PDeltaFail_neg;
                    delta_SystemFail=delta_SystemFail_neg;
                    PDeltaFail_Index=i;
                end
            end
        end
    end
    
    if Flag_PDeltaFail==1
        %uii=0.0;
        vii=0.0;
        aii=0.0;
        fii=0.0;
        kii=0.0;      
    end

    %% Save reponse quantities at end of current step
    U(i+1,1)=uii;
    V(i+1,1)=vii;
    A(i+1,1)=aii;
    F(i+1,1)=fii;
    K(i+1,1)=kii;

    %% Prepare quantitites for next step
    ui_1=ui;
    ui=uii;
    vi=vii;
    ki=kii;
    fi=fii;

    if Flag_PDeltaFail==1
        U(i+1,1)=uii;%max(abs(U(1:PDeltaFail_Index)))*1.5;
        V(i+1,1)=0.0;
        A(i+1,1)=0.0;
        F(i+1,1)=0.0;
        K(i+1,1)=0.0;      
    end
end

end