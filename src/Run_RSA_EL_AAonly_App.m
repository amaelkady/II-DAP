function [] = Run_RSA_EL_AAonly_App(app)
global MainDirectory ProjectPath ProjectName
clc
cd (ProjectPath)
load(ProjectName);
cd (MainDirectory)
clc

if exist('RESULTS')==1 && exist('RSA')==1
rmvar(ProjectPath,ProjectName, 'RESULTS', 'RSA' )
cd (MainDirectory)
end

clear RESULTS RSA

app.axes10.Position=[5 7 383 31];
app.axes10.BackgroundColor='w';
app.ProgressText.Value='RUNNING RSA ';
app.ProgressText.FontColor='y';
drawnow 

countT=1;

for TT=To:Tincr:Tend
    if TT==0; TT=0.001; end
   
for GM_No=1:nGM

    evalc(['A=','GM.GA',num2str(GM_No)]);
    evalc(['GMname=char(GM.name{',num2str(GM_No),',1})']);
    evalc(['GMpga=max(abs(GM.GA',num2str(GM_No),'))']);
    evalc(['GMdt=GM.dt(',num2str(GM_No),',1)']);

    GMpsaTi = cent_diff(A, TT, GMdt, zeta, g)/g;
    if isinf(GMpsaTi)==1; GMpsaTi=GMpga; end

    evalc(['RSA.Period(countT,1)=TT']);
    evalc(['RESULTS.AA',num2str(GM_No),'=GMpsaTi']);
    evalc(['RSA.MaxAA(countT,GM_No)=max(abs(RESULTS.AA',num2str(GM_No),'))']);



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