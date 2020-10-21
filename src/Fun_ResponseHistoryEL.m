function [Teq,U,V,A,F,K,Ag,Flag_PDeltaFail,PDeltaFail_Index] = Fun_ResponseHistoryEL (GM_No,Sa_SF, MassT,TT,CT)
global MainDirectory ProjectPath ProjectName
cd (ProjectPath)
load (ProjectName)
load (ProjectName,'gamma','beta')
cd (MainDirectory)
    
M=MassT;
T=TT;
C=CT;

nRefinePoints=GM.dt(GM_No)/dt_user;
% nRefinePoints=10;
Flag_PDeltaFail=0;
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
    Backbone.du=999;
    [kii,fii]=System_Linear (Backbone.Ke, uii,Backbone.du);
    
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

end

end