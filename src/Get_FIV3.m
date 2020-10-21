function [T, FIV]=Get_FIV3(app)
app.axes1.cla
app.axes1.FontName = 'Times';
app.axes1.FontSize = 17;

global MainDirectory ProjectPath ProjectName
clc;
cd (ProjectPath)
load(ProjectName,'GM_Option','nGM','GM','Wave_name','AddTime','g','Units','T')
cd (MainDirectory)
MaxPGA=app.axes1.UserData;
GM_No=app.GMpopup.Value;
app.text26.Tooltip='Time segment coefficient';
app.text25.Tooltip='Frequency cut-off coefficient';

T1=T;

clear T

evalc(['GMTime=','GM.Time',num2str(GM_No)]);
evalc(['Duration=','GM.Time',num2str(GM_No),'(end,1)']);
MainGMTime=Duration-AddTime;
for i=1:size(GMTime,1)
    if GMTime(i,1)>=MainGMTime
        idexMainGM=i;
        break
    end
end

evalc(['T=','GM.Time',num2str(GM_No),'(1:idexMainGM,1)']);
evalc(['GA=','GM.GA',num2str(GM_No),'(1:idexMainGM,1)']);
dt=GM.dt(GM_No,1);

beta =app.edit7.Value;
alpha=app.edit8.Value;

CutOffFreq=beta*1/T1;
TimeSegment=alpha*T1;

% Low-pass Signal Filter using the CutOffFreq
Signal=GA;
SignalFreq=1/dt;
FilterOrder=2;
[B,A] = butter(FilterOrder,CutOffFreq/(0.5*SignalFreq), 'low');
filteredSignal = transpose(filter(B, A, transpose(Signal)));

% Computing FIV history by integrating the filtered accleration signal
StepSegment=floor(TimeSegment/dt)+1;
FIV = zeros(length(filteredSignal),1);
for i=1:size(filteredSignal,1)-StepSegment
    FIV(i,1)=(sum(filteredSignal(i:i+StepSegment)))*dt;
end

% Finding the peak FIV values in both direction and computing FIV3
[MaxPosFIV,MaxPosFIVindx] = findpeaks( FIV);
[MaxNegFIV,MaxNegFIVindx] = findpeaks(-FIV);
MaxPosFIVsorted=sort(MaxPosFIV,'descend');
MaxNegFIVsorted=sort(MaxNegFIV,'descend');
FIV3=max(sum(MaxPosFIVsorted(1:3,1)),sum(abs(MaxNegFIVsorted(1:3,1))));
indxpos1=find(FIV==MaxPosFIVsorted(1,1));
indxpos2=find(FIV==MaxPosFIVsorted(2,1));
indxpos3=find(FIV==MaxPosFIVsorted(3,1));
indxneg1=find(FIV==-MaxNegFIVsorted(1,1));
indxneg2=find(FIV==-MaxNegFIVsorted(2,1));
indxneg3=find(FIV==-MaxNegFIVsorted(3,1));

app.axes1.XLabel.String='Time [sec]';
if Units==1; app.axes1.YLabel.String='V_s [m/sec]'; end
if Units==2; app.axes1.YLabel.String='V_s [in/sec]'; end
app.axes1.XLabel.FontName='Times';
app.axes1.YLabel.FontName='Times';

plot(app.axes1,T,FIV,'-r','linewidth',1.5);
scatter(app.axes1,T(indxpos1,1),FIV(indxpos1,1),'ok','markerfacecolor','b')
scatter(app.axes1,T(indxpos2,1),FIV(indxpos2,1),'ok','markerfacecolor','b')
scatter(app.axes1,T(indxpos3,1),FIV(indxpos3,1),'ok','markerfacecolor','b')
scatter(app.axes1,T(indxneg1,1),FIV(indxneg1,1),'ok','markerfacecolor','g')
scatter(app.axes1,T(indxneg2,1),FIV(indxneg2,1),'ok','markerfacecolor','g')
scatter(app.axes1,T(indxneg3,1),FIV(indxneg3,1),'ok','markerfacecolor','g')

%% Plot Origin Lines
plot(app.axes1,[0 MainGMTime*2],[0 0],'-k')
app.axes1.YLim=[-max(abs(FIV))*1.1 max(abs(FIV))*1.1];
app.axes1.XLim=[0 MainGMTime*1.1];

if Units==1
    text(app.axes1,T(indxneg3,1),max(abs(FIV)),[' FIV3=',num2str(round(FIV3*100)/100),' m/sec'], 'fontname', 'Courier', 'fontsize',12,'fontweight','bold')
else
    text(app.axes1,T(indxneg3,1),max(abs(FIV)),[' FIV3=',num2str(round(FIV3*100)/100),' in/sec'], 'fontname', 'Courier', 'fontsize',12,'fontweight','bold')
end
end