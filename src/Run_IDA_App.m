function [] = Run_IDA_App(app)
global MainDirectory ProjectPath ProjectName
clc
cd (ProjectPath)
load(ProjectName);
load(ProjectName,'gamma','beta');
cd (MainDirectory)
clc

app.axes10.Position=[5 7 383 31];
app.axes10.BackgroundColor='w';

if exist('RESULTS')==1 && exist('IDA')==1
    rmvar(ProjectPath,ProjectName, 'RESULTS', 'IDA' )
    cd (MainDirectory)
end

clear RESULTS IDA 

for GM_No=1:nGM
   
    % Set Initial Values
    count=1;
    Limit_Flag=0;
    Collapse_Flag=0;
    Sa_Current=0;
    Sa_step=Sa_incr*g;
    Sa_collapse=Alimit*g;
    if SystemModel_Option==1 || SystemModel_Option==2
        U_collapse=min([delta_SystemFail; delta_PDeltaFail; Ulimit]);
    elseif SystemModel_Option==3 || SystemModel_Option==4 || SystemModel_Option == 5 || SystemModel_Option == 6
        U_collapse=min([delta_SystemFail_pos; delta_SystemFail_neg; delta_PDeltaFail_pos; delta_PDeltaFail_neg; Ulimit]);
    end
    
    % Get Sa at T1 for the Unscaled Ground Motion
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
    
    while Collapse_Flag==0
        
        % Current Scale Factor
        Sa_SF=Sa_Current/Sa_T1;
        
        % Run Response History Analysis at current Scale Factor
        [Teq,U,V,A,F,K,Ag,Flag_PDelta_Collapse,PDeltaFail_Index] = Fun_ResponseHistory (GM_No,Sa_SF,M,T,C);
        
        U(~isfinite(U))=0;
        V(~isfinite(V))=0;
        A(~isfinite(A))=0;
        F(~isfinite(F))=0;
        K(~isfinite(K))=0;
    
        MaxU=max(abs(U));
        ResU=abs(U(end));
        MaxAA=max(abs(A+Ag));

        % Check if Collapse Criteria is Reached
        if MaxU < U_collapse  && Sa_Current < Sa_collapse          % If displacment and acceleration is less than collaspe limits
            if Sa_step > Sa_tol*g   % If Sa step is larger than the tolerance (i.e., collpase was not reached so far) 
                % Save the Current Sa and the Corrosponding Displacement
                evalc(['IDA.Sa',num2str(GM_No),'(count,1)=Sa_Current']);
                evalc(['IDA.U',num2str(GM_No),'(count,1)=MaxU']);
                evalc(['IDA.Ures',num2str(GM_No),'(count,1)=ResU']);
                evalc(['IDA.AA',num2str(GM_No),'(count,1)=MaxAA']);
                evalc(['IDA.Sa_SF',num2str(GM_No),'(count,1)=Sa_SF']);
                evalc(['IDA.U',num2str(GM_No),'_',num2str(count),'=U']);
                evalc(['IDA.F',num2str(GM_No),'_',num2str(count),'=F']);
                count=count+1;
            end
            Sa_Current = Sa_Current + Sa_step;       % Update Sa for next increment
        elseif MaxU >= U_collapse || Sa_Current >= Sa_collapse  || Limit_Flag==1
            Limit_Flag=1;
            if Sa_step > Sa_tol*g  % If collaspe displacemnt is reached and IDA step is larger than tolerance
                Sa_Current = Sa_Current - 0.5 * Sa_step; % Roll Back 1/2 Step
                Sa_step = 0.5 * Sa_step;                 % Modify Step Size to 1/2 of Previous Step 
            elseif Sa_step <= Sa_tol*g % If collaspe displacemnt is reached and IDA step is less than tolerance
                % Set Collapse Flag to 1.0 and Save the Collapse Sa and response histories at Collapse 
                Collapse_Flag=1;
                RESULTS.SaSF(GM_No,1)=Sa_SF;
                RESULTS.SaCPS(GM_No,1)=Sa_Current;
                evalc(['RESULTS.T',num2str(GM_No),'=Teq']);
                evalc(['RESULTS.U',num2str(GM_No),'=U']);
                evalc(['RESULTS.V',num2str(GM_No),'=V']);
                evalc(['RESULTS.A',num2str(GM_No),'=A']);
                evalc(['RESULTS.F',num2str(GM_No),'=F']);
                evalc(['RESULTS.K',num2str(GM_No),'=K']);
                evalc(['RESULTS.AG',num2str(GM_No),'=Ag']);
                evalc(['RESULTS.Collapse_Flag',num2str(GM_No),'=Flag_PDelta_Collapse']);
                evalc(['RESULTS.Collapse_Index',num2str(GM_No),'=PDeltaFail_Index']);
            end
        end
        
        clear AgX Teq U V A F K Ag AA AAA

        % Progress Bar
        if nGM>1
            app.axes10.Position=[5 7 (round(GM_No*100/nGM)*383/100) 31];
            app.axes10.BackgroundColor='g';
            app.ProgressText.Value=['RUNNING IDA ',num2str(round(GM_No*100/nGM)),'%'];
        else
            app.ProgressText.Value='RUNNING IDA - PLEASE WAIT...';   
            app.axes10.Position=[5 7 0.5*383 31];
            app.axes10.BackgroundColor='g';
        end
        app.ProgressText.FontColor='y';
        drawnow
        
    end
    
    
    % Save Last IDA curve data at collapse
    evalc(['IDA.Sa',num2str(GM_No),'(count,1)=Sa_Current']);
    evalc(['IDA.U',num2str(GM_No),'(count,1)=MaxU']);
    evalc(['IDA.AA',num2str(GM_No),'(count,1)=MaxAA']);
    evalc(['IDA.Sa_SF',num2str(GM_No),'(count,1)=Sa_SF']);
    evalc(['IDA.U',num2str(GM_No),'_',num2str(count),'=RESULTS.U',num2str(GM_No)]);
    evalc(['IDA.F',num2str(GM_No),'_',num2str(count),'=RESULTS.F',num2str(GM_No)]);
    % additional entry for IDA curve plot purposes
    evalc(['IDA.Sa',num2str(GM_No),'(count+1,1)=Sa_Current']);
    evalc(['IDA.U',num2str(GM_No),'(count+1,1)=MaxU*1.25']);
    evalc(['IDA.AA',num2str(GM_No),'(count+1,1)=MaxAA']);
    evalc(['IDA.Ures',num2str(GM_No),'(count,1)=ResU']);

    clear AgX Teq U V A F K Ag AA AAA Sa_Current MaxU MaxAA Sa_SF ResU
end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CALCULATE COLLAPSE FRAGILTIY
        if nGM>=3
            % Progress Bar
            app.ProgressText.Value='CALCULATING COLLAPSE FRAGILITY';
            app.ProgressText.FontColor='y';
            app.axes10.Position=[5 7 383 31];
            app.axes10.BackgroundColor='g';

            MAFc='NaN';
            Pc='NaN';
            dSa=0.005;
            clc;

            %% PLOT COLLAPSE FRAGILITY CURVE
            SaCPS=RESULTS.SaCPS/g;
            SaCPS=sort(SaCPS, 1);
            Parameters=lognfit(SaCPS);
            MedianCPS=Parameters(1);
            SigmaCPS=Parameters(2);
            if max(SaCPS)==min(SaCPS)
                MedianCPS=MedianCPS;
                SigmaCPS=0.001;
            end

            ProbabilityEmp=zeros(nGM,1);
            for i=1:nGM
                ProbabilityEmp(i,1)=i/nGM;
            end
            SAxx = (0.001:dSa:max(SaCPS)+1.0); 
            ProbabilityFit = logncdf(SAxx,MedianCPS,SigmaCPS);

            %% Calculate the mean annual frequency of collapse, MAFc
            if Hazardstatus==1

                nyears=50;

                MAFc=0.0;
                for Sa=0.001:dSa:max(SaCPS)+1.0
                    f=COEFF(1)*log(Sa)^4+COEFF(2)*log(Sa)^3+COEFF(3)*log(Sa)^2+COEFF(4)*log(Sa)^1+COEFF(5);
                    f1=COEFF(1)*log(Sa+dSa)^4+COEFF(2)*log(Sa+dSa)^3+COEFF(3)*log(Sa+dSa)^2+COEFF(4)*log(Sa+dSa)^1+COEFF(5);
                    Slope= (exp(f)-exp(f1))/dSa;
                    Probability_Fragility = logncdf(Sa,MedianCPS,SigmaCPS);
                    MAFc=MAFc+Probability_Fragility*(Slope)*dSa;
                end
                Pc=1-exp(-MAFc*nyears);
                Pc=round(Pc*10000)/10000;
            else
                nyears=50;                
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if nGM>=3
            % Save the Data
            cd (ProjectPath)
            pause(0.8)
            save(ProjectName,'RESULTS','IDA','SAxx','ProbabilityEmp','ProbabilityFit','MedianCPS','SigmaCPS','MAFc','Pc','SaCPS','MedianCPS','SigmaCPS','nyears','-append')
            pause(0.8)
            cd (MainDirectory)            
        else
            % Save the Data
            cd (ProjectPath)
            pause(0.8)
            save(ProjectName,'RESULTS','IDA','-append')
            pause(0.8)
            cd (MainDirectory)
        end
        
end