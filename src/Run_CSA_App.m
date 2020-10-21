function [] = Run_CSA_App(app)

global MainDirectory ProjectPath ProjectName
cd (ProjectPath)
load(ProjectName);
load(ProjectName,'gamma','beta');
cd (MainDirectory)
clc

nPeriods=(Tend-To)/Tincr;
nRuns=nPeriods*nGM;

clear CSA RSA IDA RESULTS

P_original=P;
H_original=H;
SC_original=Stability_Coeff;

app.ProgressText.Value='RUNNING CSA ';
app.ProgressText.FontColor='y';
drawnow 

countT=1;
    
for TT=To:Tincr:Tend

    if TT==0; TT=0.00001; end

    if RSA_Option==1
        MassT=Backbone.Ke*(TT/2/pi)^2;
    else
        MassT=M;
        [BackboneX]=UpdateRSABackbone (TT,P_original,H_original,SC_original);
    end
    omega=2*pi/TT;
    CT=zeta*2*omega*MassT;
    
    for GM_No=1:nGM
        
        count=1;

        % Set Initial Values
        Collapse_Flag=0;
        Sa_Current=0;
        Sa_step=Sa_incr*g;
        Sa_step=Sa_incr*g;
        Sa_collapse=Alimit*g;
        if SystemModel_Option==1 || SystemModel_Option==2
            U_collapse=min([delta_SystemFail; delta_PDeltaFail; Ulimit]);
        elseif SystemModel_Option==3 || SystemModel_Option==4 || SystemModel_Option == 5 || SystemModel_Option == 6
            U_collapse=min([delta_SystemFail_pos; delta_SystemFail_neg; delta_PDeltaFail_pos; delta_PDeltaFail_neg; Ulimit]);
        end
        
        % Get Sa at T1 for the Unscaled Ground Motion
        evalc([' AgX(:,1)=GM.GA',  num2str(GM_No)]);
        %Sa_T1 =  cent_diff(T, MassT, GM.dt(GM_No), zeta, AgX);
        if SaAVG_status==0
            [z1,z2,z3,AA,z4,z5,AAA,z7,z8] = Fun_ResponseHistoryEL (GM_No,1.0, MassT,TT,CT);
            Sa_T1=max(abs(AA+AAA));
            if Sa_T1==0; Sa_T1=0.001; end
        else
            Sa_PRODUCT=1;
            for Ti=SaAVG_cT1*TT:0.1:SaAVG_cT2*TT
                Mxx=Ko*(Ti/2/pi)^2;
                omega=2*pi/Ti;
                Cxx=zeta*2*omega*Mxx;
                [z1,z2,z3,AA,z4,z5,AAA,z7,z8] = Fun_ResponseHistoryEL (GM_No,1.0, Mxx,Ti,Cxx);
                Sa_Ti=max(abs(AA+AAA));
                Sa_PRODUCT=Sa_PRODUCT*Sa_Ti;
            end
            Sa_T1=(Sa_PRODUCT)^(1/11);
            if Sa_T1==0; Sa_T1=0.001; end
        end
        clear AA AAA

        while Collapse_Flag==0

            % Current Scale Factor
            Sa_SF=Sa_Current/Sa_T1;

            % Run Response History Analysis at current Scale Factor
            [Teq,U,V,A,F,K,Ag,Flag_PDelta_Collapse,PDeltaFail_Index] = Fun_ResponseHistory (GM_No,Sa_SF,MassT,TT,CT);

            U(~isfinite(U))=0;
            V(~isfinite(V))=0;
            A(~isfinite(A))=0;
            F(~isfinite(F))=0;
            K(~isfinite(K))=0;

            MaxU=max(abs(U));
            
            %disp([' GM= ',num2str(GM_No),' Sa_Current= ',num2str(Sa_Current/g),' Sa_SF= ',num2str(Sa_SF),' Sa_Step= ',num2str(Sa_step/g),' MaxU= ',num2str(MaxU)]);

            % Check if Collapse Criteria is Reached
            if MaxU < U_collapse && Sa_Current < Sa_collapse          % If displacment and acceleration is less than collaspe limits
                if Sa_step > Sa_tol*g   % If Sa step is larger than the tolerance (i.e., collpase was not reached so far) 
                    % Save the Current Sa and the Corrosponding Displacement
                    evalc(['IDA.Sa',num2str(GM_No),num2str(countT),'(count,1)=Sa_Current']);
                    evalc(['IDA.U',num2str(GM_No),num2str(countT),'(count,1)=MaxU']);
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
%                     RESULTS.SaSF(GM_No,1)=Sa_SF;
                    RESULTS.SaCPS(GM_No,countT)=Sa_Current;
%                     evalc(['RESULTS.T',num2str(GM_No),'=Teq']);
%                     evalc(['RESULTS.U',num2str(GM_No),'=U']);
%                     evalc(['RESULTS.V',num2str(GM_No),'=V']);
%                     evalc(['RESULTS.A',num2str(GM_No),'=A']);
%                     evalc(['RESULTS.F',num2str(GM_No),'=F']);
%                     evalc(['RESULTS.K',num2str(GM_No),'=K']);
%                     evalc(['RESULTS.AG',num2str(GM_No),'=Ag']);
%                     evalc(['RESULTS.Collapse_Flag',num2str(GM_No),'=Flag_PDelta_Collapse']);
%                     evalc(['RESULTS.Collapse_Index',num2str(GM_No),'=PDeltaFail_Index']);                    
                end
            end

        evalc(['RSA.Period(countT,1)=TT']);
        evalc(['RSA.U',num2str(GM_No),'_',num2str(countT),'=U']);
        evalc(['RSA.F',num2str(GM_No),'_',num2str(countT),'=F']);
        
        clear  AgX Teq U V A F K Ag AA AAA
        end
        
        % Save Last IDA curve data at collapse
        evalc(['IDA.Sa',num2str(GM_No),num2str(countT),'(count,1)=Sa_Current']);
        evalc(['IDA.U',num2str(GM_No),num2str(countT),'(count,1)=MaxU']);
        % additional entry for IDA curve plot purposes
        evalc(['IDA.Sa',num2str(GM_No),num2str(countT),'(count+1,1)=Sa_Current']);
        evalc(['IDA.U',num2str(GM_No),num2str(countT),'(count+1,1)=MaxU*1.25']);
    
        CSA.Period(countT,1)=TT;
        CSA.SaCPS(countT,GM_No)=Sa_Current;
        CSA.SaT1(countT,GM_No)=Sa_T1;
        if RSA_Option==1
            if SystemModel_Option==3 || SystemModel_Option==6
                CSA.CC(countT,GM_No)=Sa_Current/g/(Backbone.Fy_pos0/M/g);
            else
                CSA.CC(countT,GM_No)=Sa_Current/g/(Backbone.Fy_pos/M/g);                
            end
        else
            if SystemModel_Option==3 || SystemModel_Option==6
                CSA.CC(countT,GM_No)=Sa_Current/g/(BackboneX.Fy_pos0/M/g);
            else
                CSA.CC(countT,GM_No)=Sa_Current/g/(BackboneX.Fy_pos/M/g);                
            end
        end
%         evalc(['CSA.MaxU(countT,GM_No)=max(abs(RESULTS.U',num2str(GM_No),'))']);
%         evalc(['CSA.MaxV(countT,GM_No)=max(abs(RESULTS.V',num2str(GM_No),'))']);
%         evalc(['CSA.MaxA(countT,GM_No)=max(abs(RESULTS.A',num2str(GM_No),'))']);
%         evalc(['CSA.MaxF(countT,GM_No)=max(abs(RESULTS.F',num2str(GM_No),'))']);
%         evalc(['CSA.MaxK(countT,GM_No)=max(abs(RESULTS.K',num2str(GM_No),'))']);
%         evalc(['CSA.MaxAG(countT,GM_No)=max(abs(RESULTS.AG',num2str(GM_No),'))']);

    app.ProgressText.Value=['RUNNING CSA ',num2str(round((TT-To)*GM_No*100/Tincr/nRuns)),'%'];
    app.axes10.Position=[5 7 ((round((countT-1)*100/((Tend-To)/Tincr)))*383/100) 31];
    app.axes10.BackgroundColor='g';   
    app.ProgressText.FontColor='y';
    drawnow 

    end

    countT=countT+1;
 
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE COLLAPSE FRAGILTIY
    if nGM>=5
        % Progress Bar
        app.ProgressText.Value='CALCULATING COLLAPSE FRAGILITY';
        app.ProgressText.FontColor='y';
        app.axes10.Position=[5 7 383 31];
        app.axes10.BackgroundColor='g';

        MAFc='NaN';
        Pc='NaN';
        dSa=0.05;
        clc;

        %% CALCUALTE THE COLLAPSE FRAGILITY CURVE
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
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if nGM>=5
        % Save the Data
        cd (ProjectPath)
        save(ProjectName,'RSA','CSA','RESULTS','IDA','SAxx','ProbabilityEmp','ProbabilityFit','MedianCPS','SigmaCPS','MAFc','Pc','SaCPS','MedianCPS','SigmaCPS','nyears','-append')
        pause(0.8)
        cd (MainDirectory)            
    else
        % Save the Data
        cd (ProjectPath)
        save(ProjectName,'RSA','CSA','IDA','RESULTS','-append');
        pause(0.8)
        cd (MainDirectory)
    end        
    
    app.ProgressText.Value='ANALYSIS DONE';
    app.ProgressText.FontColor='g';
    app.axes10.Position=[5 7 383 31];
    app.axes10.BackgroundColor='g';

    drawnow
    
end