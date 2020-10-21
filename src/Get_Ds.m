function Get_Ds(app)
app.axes1.cla
app.axes1.FontName = 'Times';
app.axes1.FontSize = 17;

global MainDirectory ProjectPath ProjectName
clc;
cd (ProjectPath)
load(ProjectName,'GM_Option','nGM','GM','Wave_name','AddTime','g','Units')
cd (MainDirectory)
MaxPGA=app.axes1.UserData;
GM_No=app.GMpopup.Value;

clear VGx UGx Tx GAx

evalc(['Tx=','GM.Time',num2str(GM_No)]);
evalc(['GAx=','GM.GA',num2str(GM_No)]);

VGx(1,1)=0;
UGx(1,1)=0;
for i=2:size(GAx,1)
    VGx(i,1)=VGx(i-1,1)+0.5*(GAx(i,1)+GAx(i-1,1))*(Tx(i,1)-Tx(i-1,1));
    UGx(i,1)=UGx(i-1,1)+0.5*(VGx(i,1)+VGx(i-1,1))*(Tx(i,1)-Tx(i-1,1));
end

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
evalc(['GV=','VGx(1:idexMainGM,1)']);
evalc(['GU=','UGx(1:idexMainGM,1)']);
dt=GM.dt(GM_No,1);
threshold= app.edit7.Value/100;

app.axes1.XLabel.String='Time [sec]';
app.axes1.YLabel.String='Arias Intensity [%]';

app.axes1.XLabel.FontName='Times';
app.axes1.YLabel.FontName='Times';

dGA=diff(GA);
dGA_squared=dGA.^2;
dT=diff(T);
Sum=0;
for i=1:size(dT,1)
    product=pi/g*dGA_squared(i,1)*dT(i,1);
    Sum=Sum+product;
    AI(i+1,1)=Sum;
end
plot(app.axes1,T,AI*100/AI(end,1),'-r','linewidth',1.5)
app.axes1.YLim=[0 100];
app.axes1.XLim=[0 max(T)];

for i=2:size(AI,1)
    if AI(i-1,1)/AI(end,1)<=threshold && AI(i,1)/AI(end,1)>=threshold
        indexAI_start=i;
        T_at_AI_start=T(i,1);
    end
    if AI(i-1,1)/AI(end,1)<=(1-threshold) && AI(i,1)/AI(end,1)>=(1-threshold)
        indexAI_end=i;
        T_at_AI_end=T(i,1);
    end
end
%% Add label at pga
plot(app.axes1,[0 T_at_AI_start ],[threshold*100   threshold*100],'--k','linewidth',1)
plot(app.axes1,[0 T_at_AI_end],[(1-threshold)*100 (1-threshold)*100],'--k','linewidth',1)
scatter(app.axes1,[T_at_AI_start ],[(threshold)*100],'ok','Markerfacecolor','g')
scatter(app.axes1,[T_at_AI_end],[(1-threshold)*100],'ok','Markerfacecolor','g')
text(app.axes1,max(T)*0.3,50,[' Significant Duration, D_s=',num2str(abs(T_at_AI_end-T_at_AI_start)),' sec'], 'fontname', 'Courier', 'fontsize',12,'fontweight','bold')
end