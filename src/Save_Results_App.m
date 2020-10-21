function []=Save_Results_App(app)
global MainDirectory ProjectPath ProjectName

clc
[Text]=Get_Project_Summary;

cd (MainDirectory);
cd (ProjectPath);
load(ProjectName);

xx=ProjectName;
xx=xx(1:end-4);
Resultsfoldername=['Results - ' ,xx];

mkdir (Resultsfoldername);
cd(Resultsfoldername)

app.ProgressText.Value='SAVING PROJECT INFO - PLEASE WAIT';
app.ProgressText.FontColor='y';
drawnow

%% Saving Project Information
filename='Project_Summary_Data.txt';
fileX1 = fopen(filename,'wt');
for i=1:size(Text,2)
fprintf(fileX1,'%s\n',char(Text{1,i}));
end

%% RHA RESULTS
if Analysis_Option ==1
    % Save the response history data
    for GM_No=1:nGM
        zz=char(GM.name{GM_No,1});
        indexz=strfind(zz,'.');
        zz=zz(1:indexz-1);
        evalc(['filename=','''','Response History ',zz,'.txt','''']);
        fileX = fopen(filename,'wt');
        
        evalc(['SF=RESULTS.SF',num2str(GM_No)]);
        fprintf(fileX,'%s\t%f','Scale factor = ',SF);
        fprintf(fileX,'\n');
        
        if Units==1
            a='Time [sec]'; b='Ground-Acceleration [m/sec2]'; c ='Rel-Displacement [m]'; d='Rel-Velocity [m/sec]';e ='Rel-Acceleration [m/sec2]'; f='Pseudo-Velocity [m/sec]';g ='Pseudo-Acceleration [m/sec2]'; h='Abs-Acceleration [m/sec2]'; i='Force [kN]';
        else
            a='Time [sec]'; b='Ground-Acceleration [in/sec2]'; c ='Rel-Displacement [in]'; d='Rel-Velocity [in/sec]';e ='Rel-Acceleration [in/sec2]'; f='Pseudo-Velocity [in/sec]';g ='Pseudo-Acceleration [in/sec2]'; h='Abs-Acceleration [in/sec2]'; i='Force [kip]';
        end
        fprintf(fileX,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',a,b,c,d,e,f,g,h,i);
        fprintf(fileX,'\n');

        evalc(['Time=RESULTS.T',num2str(GM_No)]);
        evalc(['Disp=RESULTS.U',num2str(GM_No)]);
        evalc(['Velocity=RESULTS.V',num2str(GM_No)]);
        evalc(['Acceleration=RESULTS.A',num2str(GM_No)]);
        evalc(['PVelocity=RESULTS.PV',num2str(GM_No)]);
        evalc(['PAcceleration=RESULTS.PA',num2str(GM_No)]);
        evalc(['AbsAcceleration=RESULTS.AA',num2str(GM_No)]);
        evalc(['Force=RESULTS.F',num2str(GM_No)]);
        evalc(['GroungAcceleration=RESULTS.AG',num2str(GM_No)]);

        for i=1:size(Time,1)
          fprintf(fileX,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t', Time(i,1), GroungAcceleration(i,1), Disp(i,1), Velocity(i,1), Acceleration(i,1), PVelocity(i,1), PAcceleration(i,1),AbsAcceleration(i,1), Force(i,1));
          fprintf(fileX,'\n');
        end
        
        app.ProgressText.Value=['SAVING DYNAMIC RESULTS ',num2str(round(GM_No*100/nGM)),'%'];
        app.ProgressText.FontColor='y';
        app.axes10.Position=[5 7 (round(GM_No*100/nGM)*382/100) 31];
        app.axes10.BackgroundColor='g';
        drawnow 
    end 
end

%% IDA RESULTS
if Analysis_Option==2
    % Save the IDA curves
    for GM_No=1:nGM
        
        zz=char(GM.name{GM_No,1});
        indexz=strfind(zz,'.');
        zz=zz(1:indexz-1);
        evalc(['filename=','''','IDA DISP ',zz,'.txt','''']);
        fileX = fopen(filename,'wt');
        if Units==1
            a='Displacement [m]'; b ='Acceleration [m/sec2]';
        else
            a='Displacement [in]'; b ='Acceleration [in/sec2]';
        end
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');

        evalc(['Disp=IDA.U',num2str(GM_No)]);
        evalc(['Acceleration=IDA.Sa',num2str(GM_No)]);

        for i=1:size(Disp,1)
          fprintf(fileX,'%f\t%f\t', Disp(i,1), Acceleration(i,1));
          fprintf(fileX,'\n');
        end
        
        app.ProgressText.Value=['SAVING IDA CURVES  (IM vs Disp)',num2str(round(GM_No*100/nGM)),'%'];
        app.axes10.Position=[5 7 (round(GM_No*100/nGM)*382/100) 31];
        app.axes10.BackgroundColor='g';
        drawnow
    end
    
        % Save the IDA curves
    for GM_No=1:nGM
        
        zz=char(GM.name{GM_No,1});
        indexz=strfind(zz,'.');
        zz=zz(1:indexz-1);
        evalc(['filename=','''','IDA ResDISP ',zz,'.txt','''']);
        fileX = fopen(filename,'wt');
        if Units==1
            a='Res. Displacement [m]'; b ='Acceleration [m/sec2]';
        else
            a='Res. Displacement [in]'; b ='Acceleration [in/sec2]';
        end
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');

        evalc(['Disp=IDA.Ures',num2str(GM_No)]);
        evalc(['Acceleration=IDA.Sa',num2str(GM_No)]);

        for i=1:size(Disp,1)
          fprintf(fileX,'%f\t%f\t', Disp(i,1), Acceleration(i,1));
          fprintf(fileX,'\n');
        end
        
        app.ProgressText.Value=['SAVING IDA CURVES  (IM vs Disp)',num2str(round(GM_No*100/nGM)),'%'];
        app.axes10.Position=[5 7 (round(GM_No*100/nGM)*382/100) 31];
        app.axes10.BackgroundColor='g';
        drawnow
        
    end
    
    % Save the IDA curves
    for GM_No=1:nGM
        
        zz=char(GM.name{GM_No,1});
        indexz=strfind(zz,'.');
        zz=zz(1:indexz-1);
        evalc(['filename=','''','IDA ACC ',zz,'.txt','''']);
        fileX = fopen(filename,'wt');
        if Units==1
            if SaAVG_status==0
                a='Abs. Acceleration [m/sec2]'; b ='IM=Sa(T1) [m/sec2]';
            else
                a='Abs. Acceleration [m/sec2]'; b ='IM=Sa_avg [m/sec2]';
            end
        else
            if SaAVG_status==0
                a='Abs. Acceleration [in/sec2]'; b ='IM=Sa(T1) [m/sec2]';
            else
                a='Abs. Acceleration [in/sec2]'; b ='IM=Sa_avg [m/sec2]';
            end
        end
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');

        evalc(['AA=IDA.AA',num2str(GM_No)]);
        evalc(['Acceleration=IDA.Sa',num2str(GM_No)]);

        for i=1:size(AA,1)-2
           if AA(i,1)<2*AA(i+1,1)
          fprintf(fileX,'%f\t%f\t', AA(i,1), Acceleration(i,1));
          fprintf(fileX,'\n');
           end
        end
        
        app.ProgressText.Value=['SAVING IDA CURVES (IM vs Abs. Acc.)',num2str(round(GM_No*100/nGM)),'%'];
        app.axes10.Position=[0.0 0.076923 (round(GM_No*100/nGM)*71.5/100) 0.7];
        app.axes10.BackgroundColor='g';
        drawnow
        
    end
    
    % Save the Median, 16th and 84th percentile curves
    if nGM>5
        
        app.ProgressText.Value='SAVING MEDIAN IDA CURVE';
        drawnow
        
        [Median16_Sa, Median50_Sa, Median84_Sa, D_Vector]=Get_Median_Collapse_Fragility;
        cd (ProjectPath);
        cd(Resultsfoldername)

        evalc(['filename=','''','Median IDA Curve.txt','''']);
        fileX = fopen(filename,'wt');
        if Units==1
            if SaAVG_status==0
                a='Displacement [m]'; b ='IM=Sa(T1) [m/sec2]';
            else
                a='Displacement [m]'; b ='IM=Sa_avg [m/sec2]';
            end
        else
            if SaAVG_status==0
                a='Displacement [in]'; b ='IM=Sa(T1) [m/sec2]';
            else
                a='Displacement [in]'; b ='IM=Sa_avg [m/sec2]';
            end
        end
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');

        b='16th-Percentile'; c='Median'; d='84th-Percentile'; 
        fprintf(fileX,'%s\t%s\t%s\t',a,b,c,d);
        fprintf(fileX,'\n');

        for i=1:size(D_Vector,1)-2
            if D_Vector(i,1)<2*D_Vector(i+1,1)
                fprintf(fileX,'%f\t%f\t%f\t', D_Vector(i,1), Median16_Sa(i,1), Median50_Sa(i,1), Median84_Sa(i,1));
                fprintf(fileX,'\n');
            end
        end
        
    end
    
    % Save the Collapse Fragility curve data
    if nGM>=3
        
        app.ProgressText.Value='SAVING COLLAPSE FRAGILITY CURVE';
        drawnow
        
        % Save Empirical Fragility Curve
        filename='Collapse Fragility Curve - Empirical.txt';
        fileX = fopen(filename,'wt');
        a=['Sa(T1, ',num2str(zeta*100),'%) [g]']; b ='P|Collapse';
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');

        for i=1:size(SaCPS,1)
          fprintf(fileX,'%f\t%f\t', SaCPS(i,1), ProbabilityEmp(i,1));
          fprintf(fileX,'\n');
        end
        
        % Save Lognormal Fragiltiy Curve
        SAxx=SAxx';
        ProbabilityFit=ProbabilityFit';
        filename='Collapse Fragility Curve - Fitted Lognormal.txt';
        fileX = fopen(filename,'wt');
        a='Median Sa [g]'; b ='Dispersion';
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');
        fprintf(fileX,'%f\t%f\t', round(exp(MedianCPS)*1000)/1000, SigmaCPS);        
        fprintf(fileX,'\n');
        fprintf(fileX,'%s\t','-------------------------------');
        fprintf(fileX,'\n');
        a=['Sa(T1, ',num2str(zeta*100),'%) [g]']; b ='P|Collapse';
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');
        for i=1:size(SAxx,1)
          fprintf(fileX,'%f\t%f\t', SAxx(i,1), ProbabilityFit(i,1));
          fprintf(fileX,'\n');
        end
       
    end
    app.axes10.Position=[5 7 382 31];
    app.axes10.BackgroundColor='g';
end

%% RSA RESULTS
if Analysis_Option==3
    
    if Quick_RSA==1
        % Save the Response Spectrum curves
        for GM_No=1:nGM
            zz=char(GM.name{GM_No,1});
            indexz=strfind(zz,'.');
            zz=zz(1:indexz-1);
            evalc(['filename=','''','Response Spectra ',zz,'.txt','''']);
            fileX = fopen(filename,'wt');
            if Units==1
                a='Period [sec]'; b='Abs-Acceleration [m/sec2]';
            else
                a='Period [sec]'; b= 'Abs-Acceleration [in/sec2]';
            end
            fprintf(fileX,'%s\t%s',a,b);
            fprintf(fileX,'\n');

            Period=RSA.Period;
            AbsAcceleration=RSA.MaxAA(:,GM_No);

            for i=1:size(Period,1)
              fprintf(fileX,'%f\t%f\t', Period(i,1),AbsAcceleration(i,1));
              fprintf(fileX,'\n');
            end

            app.ProgressText.Value=['SAVING RESPONSE SPECTRA ',num2str(round(GM_No*100/nGM)),'%'];
            app.axes10.Position=[0.0 0.076923 (round(GM_No*100/nGM)*71.5/100) 0.7];
            app.axes10.BackgroundColor='g';
            drawnow
        end
    else
        % Save the Response Spectrum curves
        for GM_No=1:nGM
            zz=char(GM.name{GM_No,1});
            indexz=strfind(zz,'.');
            zz=zz(1:indexz-1);
            evalc(['filename=','''','Response Spectra ',zz,'.txt','''']);
            fileX = fopen(filename,'wt');
            if Units==1
                a='Period [sec]'; b='Rel-Displacement [m]'; c ='Rel-Velocity [m/sec]';   d='Rel-Acceleration [m/sec2]'; e ='Pseudo-Velocity [m/sec]';  f='Pseudo-Acceleration [m/sec2]';  g='Abs-Acceleration [m/sec2]';  h='Force [kN]';
            else
                a='Period [sec]'; b='Rel-Displacement [in]'; c ='Rel-Velocity [in/sec]'; d='Rel-Acceleration [in/sec2]';e ='Pseudo-Velocity [in/sec]'; f='Pseudo-Acceleration [in/sec2]'; g='Abs-Acceleration [in/sec2]'; h='Force [kip]';
            end
            fprintf(fileX,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',a,b,c,d,e,f,g,h);
            fprintf(fileX,'\n');

            Period=RSA.Period;
            Disp=RSA.MaxU(:,GM_No);
            Velocity=RSA.MaxV(:,GM_No);
            Acceleration=RSA.MaxA(:,GM_No);
            PVelocity=RSA.MaxPV(:,GM_No);
            PAcceleration=RSA.MaxPA(:,GM_No);
            AbsAcceleration=RSA.MaxAA(:,GM_No);
            Force=RSA.MaxF(:,GM_No);

            for i=1:size(Period,1)
              fprintf(fileX,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t', Period(i,1), Disp(i,1), Velocity(i,1), Acceleration(i,1), PVelocity(i,1), PAcceleration(i,1), AbsAcceleration(i,1), Force(i,1));
              fprintf(fileX,'\n');
            end

            app.ProgressText.Value=['SAVING RESPONSE SPECTRA ',num2str(round(GM_No*100/nGM)),'%'];
            app.axes10.Position=[5 7 (round(GM_No*100/nGM)*382/100) 31];
            app.axes10.BackgroundColor='g';
            drawnow
        end
    end
    
    if nGM>=3
        Yaxis=RSA.MaxAA(:,:)/g;   
        Period=RSA.Period;

        for i=1:size(Period,1)
               Sorted_Yaxis=sort(Yaxis(i,:));
               MeanY(i,1)=mean(Sorted_Yaxis);
               MedianY(i,1)=Sorted_Yaxis(1,round(nGM*0.5));
               Median16Y(i,1)=Sorted_Yaxis(1,round(nGM*0.16));
               Median84Y(i,1)=Sorted_Yaxis(1,round(nGM*0.84));
        end

        evalc(['filename=','''','Median and Mean Response Spectrum.txt','''']);
        fileX = fopen(filename,'wt');
        if Units==1
            a='Period [sec]'; b='Abs-Acceleration [m/sec2]';
        else
            a='Period [sec]'; b='Abs-Acceleration [in/sec2]';
        end
        fprintf(fileX,'%s\t%s\t',a,b);
        fprintf(fileX,'\n');

        b='16th-Percentile'; c='Median'; d='84th-Percentile'; e='Mean';
        fprintf(fileX,'%s\t%s\t%s\t%s\t',a,b,c,d,e);
        fprintf(fileX,'\n');
        
        for i=1:size(Period,1)
          fprintf(fileX,'%f\t%f\t%f\t%f\t%f\t', Period(i,1), Median16Y(i,1), MedianY(i,1), Median84Y(i,1), MeanY(i,1));
          fprintf(fileX,'\n');
        end
            
    end
end

%% CSA RESULTS
if Analysis_Option==4
    
    % Save the Response Spectrum curves
    for GM_No=1:nGM
        zz=char(GM.name{GM_No,1});
        indexz=strfind(zz,'.');
        zz=zz(1:indexz-1);
        evalc(['filename=','''','Collapse Spectra ',zz,'.txt','''']);
        fileX = fopen(filename,'wt');
        a='Period [sec]'; b='Collapse Capacity [g]'; c ='Norm. Collapse Capacity ';

        fprintf(fileX,'%s\t%s\t%s\t',a,b,c);
        fprintf(fileX,'\n');

        Period=CSA.Period;
        SaCPS=CSA.SaCPS(:,GM_No);
        CC=CSA.CC(:,GM_No);

        for i=1:size(Period,1)
          fprintf(fileX,'%f\t%f\t%f\t', Period(i,1), SaCPS(i,1), CC(i,1));
          fprintf(fileX,'\n');
        end
        
        app.ProgressText.Value=['SAVING RESPONSE SPECTRA',num2str(round(GM_No*100/nGM)),'%'];
        app.ProgressText.BackgroundColor='y';
        app.axes10.Position=[5 7 (round(GM_No*100/nGM)*382/100) 31];
        app.axes10.BackgroundColor='g';
        drawnow
    end
    
end

cd(MainDirectory)
%%
app.ProgressText.Value='DONE';
app.ProgressText.FontColor='g';
drawnow

fclose all;
end

function [Median16_Sa, Median50_Sa, Median84_Sa, D_Vector]=Get_Median_Collapse_Fragility
global MainDirectory ProjectPath ProjectName
clc;
cd (ProjectPath)
load(ProjectName,'nGM','IDA','g') 
cd (MainDirectory)

for i=1:nGM
    evalc(['Xaxes=IDA.U',num2str(i)]);
    MaxU(i,1)=max(Xaxes);
end

    counter = 1;
for D = 0.0:max(MaxU)/100:max(MaxU)

    for GM_No=1:nGM     % Loop for GMs 
        evalc(['Xaxes=IDA.U',num2str(GM_No)]);
        evalc(['Yaxes=IDA.Sa',num2str(GM_No),'/g']);
        Sa_Interp (GM_No,1)=interp1(Xaxes,Yaxes,D);
    end

    Sa_Interp=sort(Sa_Interp); 
    Median16_Sa(counter,1)=(Sa_Interp(ceil(0.16*nGM),1)+Sa_Interp(floor(0.16*nGM),1))/2;  
    Median50_Sa(counter,1)=median(Sa_Interp);
    Median84_Sa(counter,1)=(Sa_Interp(ceil(0.84*nGM),1)+Sa_Interp(floor(0.84*nGM),1))/2;  
    D_Vector (counter,1)=D;

    counter = counter + 1;
end

    Median16_Sa(counter,1)=Median16_Sa(end,1);  
    Median50_Sa(counter,1)=Median50_Sa(end,1);
    Median84_Sa(counter,1)=Median84_Sa(end,1);  
    Median16_Sa(~isfinite(Median16_Sa))=max(Median16_Sa);
    Median50_Sa(~isfinite(Median50_Sa))=max(Median50_Sa);
    Median84_Sa(~isfinite(Median84_Sa))=max(Median84_Sa);

    D_Vector (counter,1)=max(MaxU);
    
end
