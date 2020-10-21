function [Text]=Get_Project_Summary
global MainDirectory ProjectPath ProjectName

clc
cd (MainDirectory);
cd (ProjectPath);
%load(ProjectName,'Units','SDoFstatus','PDeltastatus','SystemModelstatus','SystemModel_Option','SymmetricBackbone','GMstatus','AnalysisTypestatus','Analysis_Option','Hazardstatus','Hazard_Option','Integratorstatus','Integrator_Option','Ko','T','M','C','zeta','Backbone','gamma','beta','dt_user')
load(ProjectName)
load(ProjectName,'gamma','beta')
cd (MainDirectory);

if Units==1
    massUnit=' kN.sec2/m';
    dispUnit=' m';
    forceUnit=' kN';
    KUnit=' kN/m';
else
    massUnit=' kip.sec2/in';    
    dispUnit=' in';
    forceUnit=' kip';
    KUnit=' kip/in';
end

i=1;
Text{i}='******************************************************';i=i+1;
Text{i}='******************************************************';i=i+1;
Text{i}='                    PROJECT SUMMARY                   ';i=i+1;
Text{i}='******************************************************';i=i+1;
Text{i}='******************************************************';i=i+1;
Text{i}='';i=i+1;

if SDoFstatus~=1
    Text{i}='NO PROJECT DATA ARE DEFINED YET';i=i+1;    
end

if SDoFstatus==1
    Text{i}='SDoF PROPERTIES';i=i+1;
    Text{i}='------------------------------------------------------';i=i+1;
    Text{i}=['Period: ', num2str(T), ' sec'];i=i+1;
    Text{i}=['Mass: ', num2str(M), massUnit];i=i+1;
    Text{i}=['Stiffness: ', num2str(round(Ko*1000)/1000), KUnit];i=i+1;
    Text{i}=['Damping Coefficient: ', num2str(zeta*100), '%'];i=i+1;
    if PDeltastatus==1
        Text{i}=['P-Delta Inclusion: ', 'YES'];i=i+1;
        Text{i}=['Stability Coefficient: ', num2str(Stability_Coeff)];i=i+1;
    else
        Text{i}=['P-Delta Inclusion: ', 'NO'];i=i+1;
    end
    Text{i}='------------------------------------------------------';i=i+1;
end

Text{i}='';i=i+1;

if SystemModelstatus==1
    Text{i}='SYSTEM TYPE';i=i+1;
    Text{i}='------------------------------------------------------';i=i+1;
    if SystemModel_Option==1
        Text{i}=['Linear System Model'];i=i+1;
        Text{i}=['Stiffness: ', num2str(round(Ko*1000)/1000), KUnit];i=i+1;
    elseif SystemModel_Option==2
        Text{i}=['Bi-Linear System Model'];i=i+1;
        Text{i}=['Fy: ', num2str(round(BackboneNoPD.Fy*1000)/1000), forceUnit];i=i+1;
        Text{i}=['dy: ', num2str(round(BackboneNoPD.dy*1000)/1000), dispUnit];i=i+1;
        Text{i}=['Ke: ', num2str(round(Ko*1000)/1000), KUnit];i=i+1;
        Text{i}=['Kp: ', num2str(BackboneNoPD.Kp), KUnit];i=i+1;
        Text{i}=['Ultimate disp.: ', num2str(delta_SystemFail), dispUnit];i=i+1;
    elseif SystemModel_Option==3
        Text{i}=['IMK-Bilinear Model'];i=i+1;
        if SymmetricBackbone==0
            Text{i}=['Assymmetric Backbone'];i=i+1;
            Text{i}=['Ke:   ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy+:  ', num2str(round(BackboneNoPD.Uy_pos0*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp+:  ', num2str(round(BackboneNoPD.Up_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc+: ', num2str(round(BackboneNoPD.Upc_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['du+:  ', num2str(round(BackboneNoPD.Uu_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy+:  ', num2str(round(BackboneNoPD.Fy_pos0*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos0/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy+: ', num2str(BackboneNoPD.FmaxFy_pos0)];i=i+1;
            Text{i}=['Fres/Fy+: ', num2str(BackboneNoPD.FresFy_pos0)];i=i+1;
            Text{i}=['dy-:  ', num2str(round(BackboneNoPD.Uy_neg0*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp-:  ', num2str(round(BackboneNoPD.Up_neg0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_neg0/BackboneNoPD.Uy_neg0*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc-: ', num2str(round(BackboneNoPD.Upc_neg0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_neg0/BackboneNoPD.Uy_neg0*100)/100),'*dy)'];i=i+1;
            Text{i}=['du-:  ', num2str(round(BackboneNoPD.Uu_neg0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_neg0/BackboneNoPD.Uy_neg0*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy-:  ', num2str(round(BackboneNoPD.Fy_neg0*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_neg0/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy-: ', num2str(BackboneNoPD.FmaxFy_neg0)];i=i+1;
            Text{i}=['Fres/Fy-: ', num2str(BackboneNoPD.FresFy_neg0)];i=i+1;
        else
            Text{i}=['Symmetric Backbone'];i=i+1;
            Text{i}=['Ke:  ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy:  ', num2str(round(BackboneNoPD.Uy_pos0*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp:  ', num2str(round(BackboneNoPD.Up_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc: ', num2str(round(BackboneNoPD.Upc_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['du:  ', num2str(round(BackboneNoPD.Uu_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy:  ', num2str(round(BackboneNoPD.Fy_pos0*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos0/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy: ', num2str(BackboneNoPD.FmaxFy_pos0)];i=i+1;
            Text{i}=['Fres/Fy: ', num2str(BackboneNoPD.FresFy_pos0)];i=i+1;
        end
    elseif SystemModel_Option==4
            Text{i}=['IMK-Pinching Model'];i=i+1;            
        if SymmetricBackbone==0
            Text{i}=['Assymmetric Backbone'];i=i+1;
            Text{i}=['Ke:   ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy+:  ', num2str(round(BackboneNoPD.Uy_pos*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp+:  ', num2str(round(BackboneNoPD.Up_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc+: ', num2str(round(BackboneNoPD.Upc_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['du+:  ', num2str(round(BackboneNoPD.Uu_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy+:  ', num2str(round(BackboneNoPD.Fy_pos*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy+: ', num2str(BackboneNoPD.FmaxFy_pos)];i=i+1;
            Text{i}=['Fres/Fy+: ', num2str(BackboneNoPD.FresFy_pos)];i=i+1;
            Text{i}=['dy-:  ', num2str(round(BackboneNoPD.Uy_neg*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp-:  ', num2str(round(BackboneNoPD.Up_neg*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_neg/BackboneNoPD.Uy_neg*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc-: ', num2str(round(BackboneNoPD.Upc_neg*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_neg/BackboneNoPD.Uy_neg*100)/100),'*dy)'];i=i+1;
            Text{i}=['du-:  ', num2str(round(BackboneNoPD.Uu_neg*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_neg/BackboneNoPD.Uy_neg*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy-:  ', num2str(round(BackboneNoPD.Fy_neg*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_neg/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy-: ', num2str(BackboneNoPD.FmaxFy_neg)];i=i+1;
            Text{i}=['Fres/Fy-: ', num2str(BackboneNoPD.FresFy_neg)];i=i+1;
            Text{i}=['kFORCE: ', num2str(BackboneNoPD.kappaF)];i=i+1;
            Text{i}=['kDISP:  ', num2str(BackboneNoPD.kappaD)];i=i+1;
        else
            Text{i}=['Symmetric Backbone'];i=i+1;
            Text{i}=['Ke:  ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy:  ', num2str(round(BackboneNoPD.Uy_pos*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp:  ', num2str(round(BackboneNoPD.Up_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc: ', num2str(round(BackboneNoPD.Upc_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['du:  ', num2str(round(BackboneNoPD.Uu_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy:  ', num2str(round(BackboneNoPD.Fy_pos*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy: ', num2str(BackboneNoPD.FmaxFy_pos)];i=i+1;
            Text{i}=['Fres/Fy: ', num2str(BackboneNoPD.FresFy_pos)];i=i+1;
            Text{i}=['kFORCE: ', num2str(BackboneNoPD.kappaF)];i=i+1;
            Text{i}=['kDISP:  ', num2str(BackboneNoPD.kappaD)];i=i+1;
        end
    elseif SystemModel_Option==5
            Text{i}=['IMK-Peak Oriented Model'];i=i+1;
        if SymmetricBackbone==0
            Text{i}=['Assymmetric Backbone'];i=i+1;
            Text{i}=['Ke:   ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy+:  ', num2str(round(BackboneNoPD.Uy_pos*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp+:  ', num2str(round(BackboneNoPD.Up_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc+: ', num2str(round(BackboneNoPD.Upc_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['du+:  ', num2str(round(BackboneNoPD.Uu_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy+:  ', num2str(round(BackboneNoPD.Fy_pos*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy+: ', num2str(BackboneNoPD.FmaxFy_pos)];i=i+1;
            Text{i}=['Fres/Fy+: ', num2str(BackboneNoPD.FresFy_pos)];i=i+1;
            Text{i}=['dy-:  ', num2str(round(BackboneNoPD.Uy_neg*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp-:  ', num2str(round(BackboneNoPD.Up_neg*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_neg/BackboneNoPD.Uy_neg*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc-: ', num2str(round(BackboneNoPD.Upc_neg*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_neg/BackboneNoPD.Uy_neg*100)/100),'*dy)'];i=i+1;
            Text{i}=['du-:  ', num2str(round(BackboneNoPD.Uu_neg*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_neg/BackboneNoPD.Uy_neg*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy-:  ', num2str(round(BackboneNoPD.Fy_neg*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_neg/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy-: ', num2str(BackboneNoPD.FmaxFy_neg)];i=i+1;
            Text{i}=['Fres/Fy-: ', num2str(BackboneNoPD.FresFy_neg)];i=i+1;
            Text{i}=['kFORCE: ', num2str(BackboneNoPD.kappaF)];i=i+1;
            Text{i}=['kDISP:  ', num2str(BackboneNoPD.kappaD)];i=i+1;
        else
            Text{i}=['Symmetric Backbone'];i=i+1;
            Text{i}=['Ke:  ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy:  ', num2str(round(BackboneNoPD.Uy_pos*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp:  ', num2str(round(BackboneNoPD.Up_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc: ', num2str(round(BackboneNoPD.Upc_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['du:  ', num2str(round(BackboneNoPD.Uu_pos*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos/BackboneNoPD.Uy_pos*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy:  ', num2str(round(BackboneNoPD.Fy_pos*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy: ', num2str(BackboneNoPD.FmaxFy_pos)];i=i+1;
            Text{i}=['Fres/Fy: ', num2str(BackboneNoPD.FresFy_pos)];i=i+1;
            Text{i}=['kFORCE: ', num2str(BackboneNoPD.kappaF)];i=i+1;
            Text{i}=['kDISP:  ', num2str(BackboneNoPD.kappaD)];i=i+1;
        end
    elseif SystemModel_Option==6
            Text{i}=['Flag-Shaped Model'];i=i+1;
        if SymmetricBackbone==0
            Text{i}=['Assymmetric Backbone'];i=i+1;
            Text{i}=['Ke:   ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy+:  ', num2str(round(BackboneNoPD.Uy_pos0*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp+:  ', num2str(round(BackboneNoPD.Up_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc+: ', num2str(round(BackboneNoPD.Upc_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['du+:  ', num2str(round(BackboneNoPD.Uu_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy+:  ', num2str(round(BackboneNoPD.Fy_pos0*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos0/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy+: ', num2str(BackboneNoPD.FmaxFy_pos0)];i=i+1;
            Text{i}=['Flim/Fy+: ', num2str(BackboneNoPD.FllFy_pos0)];i=i+1;
            Text{i}=['Kul/Ke+: ', num2str(BackboneNoPD.KulKe_npos0)];i=i+1;
            Text{i}=['dy-:  ', num2str(round(BackboneNoPD.Uy_neg0*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp-:  ', num2str(round(BackboneNoPD.Up_neg0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_neg0/BackboneNoPD.Uy_neg0*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc-: ', num2str(round(BackboneNoPD.Upc_neg0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_neg0/BackboneNoPD.Uy_neg0*100)/100),'*dy)'];i=i+1;
            Text{i}=['du-:  ', num2str(round(BackboneNoPD.Uu_neg0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_neg0/BackboneNoPD.Uy_neg0*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy-:  ', num2str(round(BackboneNoPD.Fy_neg0*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_neg0/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy-: ', num2str(BackboneNoPD.FmaxFy_neg0)];i=i+1;
            Text{i}=['Flim/Fy-: ', num2str(BackboneNoPD.FllFy_neg0)];i=i+1;
            Text{i}=['Kul/Ke-: ', num2str(BackboneNoPD.KulKe_neg0)];i=i+1;
        else
            Text{i}=['Symmetric Backbone'];i=i+1;
            Text{i}=['Ke:  ', num2str(round(BackboneNoPD.Ke*1000)/1000), KUnit];i=i+1;
            Text{i}=['dy:  ', num2str(round(BackboneNoPD.Uy_pos0*1000)/1000), dispUnit];i=i+1;
            Text{i}=['dp:  ', num2str(round(BackboneNoPD.Up_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Up_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['dpc: ', num2str(round(BackboneNoPD.Upc_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Upc_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['du:  ', num2str(round(BackboneNoPD.Uu_pos0*1000)/1000), dispUnit,' (', num2str(round(BackboneNoPD.Uu_pos0/BackboneNoPD.Uy_pos0*100)/100),'*dy)'];i=i+1;
            Text{i}=['Fy:  ', num2str(round(BackboneNoPD.Fy_pos0*1000)/1000), forceUnit,' (', num2str(round(BackboneNoPD.Fy_pos0/M/g*100)/100),'*W)'];i=i+1;
            Text{i}=['Fmax/Fy: ', num2str(BackboneNoPD.FmaxFy_pos0)];i=i+1;
            Text{i}=['Flim/Fy: ', num2str(BackboneNoPD.FllFy_pos0)];i=i+1;
            Text{i}=['Kul/Ke: ', num2str(BackboneNoPD.KulKe_pos0)];i=i+1;
        end
    end
    
    Text{i}='------------------------------------------------------';i=i+1;
end

Text{i}='';i=i+1;

if GMstatus==1
    Text{i}=['GROUND MOTION'];i=i+1;
    Text{i}='------------------------------------------------------';i=i+1;
    if GM_Option==1
        Text{i}=['Ground Motion Type: Single Ground Motion Record'];i=i+1;
        Text{i}=['Peak Ground Acceleration: ', num2str(round(GM.pga*1000)/1000), ' g'];i=i+1;    
        Text{i}=['GM Time Step: ', num2str(GM.dt), ' sec'];i=i+1;    
        if AddTime~=0
            Text{i}=['GM Duration: ', num2str(GM.duration-AddTime), ' sec + ',num2str(AddTime),' sec free vibration' ];i=i+1;    
        else
            Text{i}=['GM Duration: ', num2str(GM.duration), ' sec'];i=i+1;    
        end
        Text{i}=['Scale Factor: ', num2str(SF)];i=i+1;    
    elseif GM_Option==2
        Text{i}=['Ground Motion Type: Multiple Ground Motion Records'];i=i+1;
        Text{i}=['Number of Ground Motions: ', num2str(nGM)];i=i+1; 
        Text{i}=['Scale Factor: ', num2str(SF)];i=i+1;         
        Text{i}=['Free Vibration Time: ', num2str(AddTime)];i=i+1;         
    elseif GM_Option==3
        Text{i}=['Ground Motion Type: ',Wave_name];i=i+1;
        Text{i}=['Wave Amplitude: ', num2str(round(GM.pga*1000)/1000), ' g'];i=i+1;    
        Text{i}=['Wave Period: ', num2str(Twave), ' sec'];i=i+1;   
        if AddTime~=0
            Text{i}=['Wave Duration: ', num2str(GM.Time1(end,1)-AddTime), ' sec + ',num2str(AddTime),' sec free vibration' ];i=i+1;    
        else
            Text{i}=['Wave Duration: ', num2str(GM.Time1(end,1)), ' sec'];i=i+1;    
        end
    end
    Text{i}='------------------------------------------------------';i=i+1;
end

Text{i}='';i=i+1;

if AnalysisTypestatus==1
    Text{i}=['ANALYSIS TYPE'];i=i+1;
    Text{i}='------------------------------------------------------';i=i+1;
    if Analysis_Option==1
        Text{i}='Analysis Type: Response History Analysis';i=i+1;
        if RHA_Option==1
            Text{i}='Scaling Option: Constant Scale Factor';i=i+1;
            Text{i}=['Scale Factor= ',num2str(SF)];i=i+1;
        else
            if SaAVG_status==0
                Text{i}='Scaling Option: Target Spectral Acceleration at T1';i=i+1;
                Text{i}=['Target SA= ',num2str(TargetSA),'g'];i=i+1;   
            else
                Text{i}='Scaling Option: Target Spectral Acceleration at Sa,avg';i=i+1;
                Text{i}=['Target SA= ',num2str(TargetSA),'g'];i=i+1;       
                Text{i}=['Sa,avg period range= ',num2str(SaAVG_T1/T),'T to ',num2str(SaAVG_T2/T),'T'];i=i+1;       
            end        
        end
    elseif Analysis_Option==2
        Text{i}='Analysis Type: Incremental Dynamic Analysis';i=i+1;
        Text{i}=['IDA Step Increment: ',num2str(Sa_incr),' g'];i=i+1;
        Text{i}=['IDA Step Tolerance: ',num2str(Sa_tol),' g'];i=i+1;
        Text{i}=['Collapse Disp. Limit: ',num2str(Ulimit),dispUnit];i=i+1;
        Text{i}=['Collapse Acc. Limit: ',num2str(Alimit),' g'];i=i+1;
        if SaAVG_status==0
            Text{i}='Scaling Option: Spectral Acceleration at T1';i=i+1; 
        else
            Text{i}='Scaling Option: Spectral Acceleration at Sa,avg';i=i+1;
            Text{i}=['Sa,avg period range= ',num2str(SaAVG_T1/T),'T to ',num2str(SaAVG_T2/T),'T'];i=i+1;       
        end       
    elseif Analysis_Option==3
        Text{i}='Analysis Type: Response Spectrum Analysis';i=i+1;
        Text{i}=['Start Period: ',num2str(To),' sec'];i=i+1;
        Text{i}=['Increment Period: ',num2str(Tincr),' sec'];i=i+1;
        Text{i}=['End Period: ',num2str(Tend),' sec'];i=i+1;
    elseif Analysis_Option==4
        Text{i}='Analysis Type: Collapse Spectrum Analysis';i=i+1;
        Text{i}=['Start Period: ',num2str(To),' sec'];i=i+1;
        Text{i}=['Increment Period: ',num2str(Tincr),' sec'];i=i+1;
        Text{i}=['End Period: ',num2str(Tend),' sec'];i=i+1;
        Text{i}=['IDA Step Increment: ',num2str(Sa_incr),' g'];i=i+1;
        Text{i}=['Collapse Point Tolerance: ',num2str(Sa_tol),' g'];i=i+1;
        Text{i}=['Collapse Disp. Limit: ',num2str(Ulimit),dispUnit];i=i+1;
        Text{i}=['Collapse Acc. Limit: ',num2str(Alimit),' g'];i=i+1;
    end
    Text{i}='------------------------------------------------------';i=i+1;
end

Text{i}='';i=i+1;

if Hazardstatus==1 && Hazard_Option==1
    Text{i}=['SEISMIC HAZARD'];i=i+1;
    Text{i}='------------------------------------------------------';i=i+1;
    Text{i}='USGS Hazard Maps';i=i+1;
    Text{i}=['Longitude: ',num2str(Longitude), ' degrees'];i=i+1;
    Text{i}=['Latitude: ',num2str(Latitude), ' degrees'];i=i+1;
    Text{i}=['Period: ',num2str(T_hazard),' sec'];i=i+1;
    Text{i}=['Soil Shear Wave Velocity: ',num2str(Vs30),' m/sec'];i=i+1;
    Text{i}='------------------------------------------------------';i=i+1;
    Text{i}='';i=i+1;
end


if Integratorstatus==1
    Text{i}='INTEGRATOR TYPE';i=i+1;
    Text{i}='------------------------------------------------------';i=i+1;
    if Integrator_Option==1
        Text{i}=['Integrator: Central Difference Method'];i=i+1;
        Text{i}=['Integration Time Step: ',num2str(dt_user), ' sec'];i=i+1;
    elseif Integrator_Option==2
        Text{i}=['Integrator: Newmark'];i=i+1;
        Text{i}=['gamma: ',num2str(gamma)];i=i+1;
        Text{i}=['beta: ',num2str(beta)];i=i+1;
        Text{i}=['Integration Time Step: ',num2str(dt_user), ' sec'];i=i+1;
    elseif Integrator_Option==3
        Text{i}=['Integrator Type: Recurrrence Formulas'];i=i+1;
        Text{i}=['Integration Time Step: ',num2str(dt_user), ' sec'];i=i+1;
    end
    Text{i}=['Time Step: ',num2str(dt_user),' sec'];
end

Text{i}='******************************************************';i=i+1;
Text{i}='******************************************************';i=i+1;
