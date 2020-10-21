function [] = Run_RSA_App(app)
global MainDirectory ProjectPath ProjectName
clc
cd (ProjectPath)
load(ProjectName);
load(ProjectName,'gamma','beta');
cd (MainDirectory)
clc

if exist('RESULTS')==1 && exist('RSA')==1
rmvar(ProjectPath,ProjectName, 'RESULTS', 'RSA' )
cd (MainDirectory)
end

clear RESULTS RSA

P_original=P;
H_original=H;
SC_original=Stability_Coeff;

app.axes10.Position=[5 7 383 31];
app.axes10.BackgroundColor='w';
app.ProgressText.Value='RUNNING RSA ';
app.ProgressText.FontColor='y';
drawnow 
    
    countT=1;
    
for TT=To:Tincr:Tend
    if TT==0; TT=0.001; end
    
    if RSA_Option==1
        MassT=Backbone.Ke*(TT/2/pi)^2;
    else
        MassT=M;
        [BackboneX]=UpdateRSABackbone (TT,P_original,H_original,SC_original);
    end
    omega=2*pi/TT;
    CT=zeta*2*omega*MassT;

    for GM_No=1:nGM

        [Teq,U,V,A,F,K,Ag,Flag_PDelta_Collapse,PDeltaFail_Index] = Fun_ResponseHistory (GM_No,1.0, MassT,TT,CT);

        U(~isfinite(U))=0;
        V(~isfinite(V))=0;
        A(~isfinite(A))=0;
        F(~isfinite(F))=0;
        K(~isfinite(K))=0;

        evalc(['RESULTS.T',num2str(GM_No),'=Teq']);
        evalc(['RESULTS.U',num2str(GM_No),'=U']);
        evalc(['RESULTS.V',num2str(GM_No),'=V']);
        evalc(['RESULTS.A',num2str(GM_No),'=A']);
        evalc(['RESULTS.F',num2str(GM_No),'=F']);
        evalc(['RESULTS.K',num2str(GM_No),'=K']);
        evalc(['RESULTS.AG',num2str(GM_No),'=Ag']);
        evalc(['RESULTS.AA',num2str(GM_No),'=A+Ag']);
        evalc(['RESULTS.PV',num2str(GM_No),'=U*omega']);
        evalc(['RESULTS.PA',num2str(GM_No),'=U*omega^2']);
        evalc(['RESULTS.Collapse_Flag',num2str(GM_No),'=Flag_PDelta_Collapse']);

        evalc(['RSA.Period(countT,1)=TT']);
        evalc(['RSA.MaxU(countT,GM_No)=max(abs(RESULTS.U',num2str(GM_No),'))']);
        evalc(['RSA.MaxV(countT,GM_No)=max(abs(RESULTS.V',num2str(GM_No),'))']);
        evalc(['RSA.MaxA(countT,GM_No)=max(abs(RESULTS.A',num2str(GM_No),'))']);
        evalc(['RSA.MaxF(countT,GM_No)=max(abs(RESULTS.F',num2str(GM_No),'))']);
        evalc(['RSA.MaxK(countT,GM_No)=max(abs(RESULTS.K',num2str(GM_No),'))']);
        evalc(['RSA.MaxAG(countT,GM_No)=max(abs(RESULTS.AG',num2str(GM_No),'))']);
        evalc(['RSA.MaxAA(countT,GM_No)=max(abs(RESULTS.AA',num2str(GM_No),'))']);
        evalc(['RSA.MaxPA(countT,GM_No)=max(abs(RESULTS.PA',num2str(GM_No),'))']);
        evalc(['RSA.MaxPV(countT,GM_No)=max(abs(RESULTS.PV',num2str(GM_No),'))']);

        evalc(['RSA.U',num2str(GM_No),'_',num2str(countT),'=U']);
        evalc(['RSA.F',num2str(GM_No),'_',num2str(countT),'=F']);
        
        clear  Teq U V A F K Ag

    end

    app.axes10.Position=[5 7 ((round((countT-1)*100/((Tend-To)/Tincr)))*383/100) 31];
    app.axes10.BackgroundColor='g';
    app.ProgressText.Value=['RUNNING RSA ',num2str(round((countT-1)*100/((Tend-To)/Tincr))),'%'];
    drawnow 
    
    countT=countT+1;
 
end

    cd (ProjectPath)
    save(ProjectName,'RESULTS','RSA','-append');
    pause(0.5)
    cd (MainDirectory)
    
    app.ProgressText.Value='ANALYSIS DONE';
    app.ProgressText.FontColor='g';
    app.axes10.Position=[5 7 383 31];
    app.axes10.BackgroundColor='g';
    drawnow
    
end