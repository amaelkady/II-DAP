function [] = Run_RHA_App(app)
global MainDirectory ProjectPath ProjectName
clc
cd (ProjectPath)
load (ProjectName,'RHA_Option','RESULTS','M','T','C','GM','nGM','SF','TargetSA','g','SaAVG_status','SaAVG_T1','SaAVG_T2','Ko','zeta','SaAVG_sample')
cd (MainDirectory)

if exist('RESULTS')==1
rmvar(ProjectPath,ProjectName, 'RESULTS')
cd (MainDirectory)
end

clear RESULTS
omega=2*pi/T;

for GM_No=1:nGM

    if nGM>1
        app.axes10.Position=[5 7 (round(GM_No*100/nGM)*383/100) 31];
        app.axes10.BackgroundColor='g';
        app.ProgressText.Value=['RUNNING RHA ',num2str(round(GM_No*100/nGM)),'%'];
    else
        app.ProgressText.Value='RUNNING RHA - PLEASE WAIT...'; 
        app.axes10.Position=[5 7 0.5*383 31];
        app.axes10.BackgroundColor='g';
    end
    app.ProgressText.FontColor='y';
    drawnow
    
    if RHA_Option==1
       SF=SF; 
    else
         if SaAVG_status==0
            [z1,z2,z3,AA,z4,z5,AAA,z7,z8] = Fun_ResponseHistoryEL (GM_No,1.0, M,T,C);
            Sa_T1=max(abs(AA+AAA));
            if Sa_T1==0; Sa_T1=0.001; end
        else
            Sa_PRODUCT=1;
            nsample=0;
            sampleStep=(SaAVG_T2-SaAVG_T1)/SaAVG_sample;
            for Ti=SaAVG_T1:sampleStep:SaAVG_T2
                nsample=nsample+1;
            end
            for Ti=SaAVG_T1:sampleStep:SaAVG_T2
                Mxx=Ko*(Ti/2/pi)^2;
                omega=2*pi/Ti;
                Cxx=zeta*2*omega*Mxx;
                [z1,z2,z3,AA,z4,z5,AAA,z7,z8] = Fun_ResponseHistoryEL (GM_No,1.0, Mxx,Ti,Cxx);
                Sa_Ti=max(abs(AA+AAA));
                Sa_PRODUCT=Sa_PRODUCT*Sa_Ti^(1/nsample);
            end
            Sa_T1=(Sa_PRODUCT);
            if Sa_T1==0; Sa_T1=0.001; end
        end
        
        SF=TargetSA/(Sa_T1/g); 
        clear z1 z2 z3 AA z4 z5 AAA z7 z8
    end
    
    [Teq,U,V,A,F,K,Ag,Flag_PDelta_Collapse,PDeltaFail_Index] = Fun_ResponseHistory (GM_No,SF,M,T,C);
    
    U(~isfinite(U))=0;
    V(~isfinite(V))=0;
    A(~isfinite(A))=0;
    F(~isfinite(F))=0;
    K(~isfinite(K))=0;
      
    evalc(['RESULTS.SF',num2str(GM_No),'=SF']);
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
    evalc(['RESULTS.Collapse_Index',num2str(GM_No),'=PDeltaFail_Index']);
    
    clear Teq U V A F K Ag

end

    cd (ProjectPath)
    save(ProjectName,'RESULTS','-append');
    pause(0.5)
    cd (MainDirectory)

    app.ProgressText.Value='ANALYSIS DONE';
    app.ProgressText.FontColor='g';
    app.axes10.Position=[5 7 383 31];
    app.axes10.BackgroundColor='g';

    drawnow
    
end