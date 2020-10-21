function []=Show_Project_Summary
clc

[Text]=Get_Project_Summary;

Project_Summary
ayy=findall(groot, 'Name', 'Project Summary'); %myguide1 is the file name for guide1
set(findobj(ayy(1),'Tag','SummaryText'),'String',Text);
drawnow
