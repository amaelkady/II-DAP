function [HazardCurveData]= Get_HazardCurve_USGSMaps (HazardFolderPath,Latitude,Longitude,T,Vs30)

IMstart=0.001;
IMincr=0.002;
IMend=5.0;

cd (HazardFolderPath)

% Discrete periods given by the USGS
Period=[0.0 10 20 30 50 75 100 200 300 400 500];

% Get file indeces of preceeding and succedding USGS periods 
for i=1:length(Period)-1
    if T*100>=Period(i) && T*100<Period(i+1)
       indexPeriod=i; 
       break
    end
end

% Go through the two USGS hazard data files
for i=indexPeriod:indexPeriod+1
    if  Period(i)>=100
        evalc(['load(','''','2008.WUS.',num2str(Period(i)),'.',num2str(Vs30),'.mat','''',')']);
    else
        evalc(['load(','''','2008.WUS.0',num2str(Period(i)),'.',num2str(Vs30),'.mat','''',')']);
    end

    % Get the SA and the AFE values
    for j=1:size(Untitled,1)
        if isnan(Untitled(j,3))
        else
            indexSAend=j-1;
            break;
        end
    end
        
    SArange=Untitled(1:indexSAend,1); % Sa values in the hazard files
    AFE = Untitled(indexSAend+1:end,:); % all the AFE values at all coordinates

    % Loop over all the AFE values
    for j=1:size(AFE,1)
        % Find the row that corrospond to the building coordinates
        if abs(Latitude-AFE(j,1)) < 0.05 && abs(Longitude-AFE(j,2)) <0.05
            % select the preceeding or succeding row number based on the
            % difference in Long & Lat
            if abs(Longitude-AFE(j,2))>abs(Longitude-AFE(j+1,2))
                row=j+1;
            else
                row=j;                
            end
            break;
        end   
    end
    
    if i==indexPeriod
        % Save the discrete AFE data for the building coordinates
        SA1=SArange;
        AFE1=AFE(row,3:end); AFE1 =AFE1';
        T1=Period(i);
        Hazard_Fit_COEFF1=polyfit(log(SA1),log(AFE1),4);
    else
        % Save the discrete SA and AFE data for the building coordinates
        SA2=SArange;    
        AFE2=AFE(row,3:end); AFE2 =AFE2';
        T2=Period(i);
        Hazard_Fit_COEFF2=polyfit(log(SA2),log(AFE2),4);
    end
    
    clear SA AFE Untitled;
end

% Get AFE data at all IM points
SArange=IMstart:IMincr:IMend;
SArange=SArange';
for i=1:size(SArange,1)
    AFE1(i,1)=Hazard_Fit_COEFF1(1,1)*log(SArange(i,1))^4+Hazard_Fit_COEFF1(1,2)*log(SArange(i,1))^3+Hazard_Fit_COEFF1(1,3)*log(SArange(i,1))^2+Hazard_Fit_COEFF1(1,4)*log(SArange(i,1))^1+Hazard_Fit_COEFF1(1,5);
    AFE2(i,1)=Hazard_Fit_COEFF2(1,1)*log(SArange(i,1))^4+Hazard_Fit_COEFF2(1,2)*log(SArange(i,1))^3+Hazard_Fit_COEFF2(1,3)*log(SArange(i,1))^2+Hazard_Fit_COEFF2(1,4)*log(SArange(i,1))^1+Hazard_Fit_COEFF2(1,5);
end

% Interpolate to get values at building T1
for i=1:size(AFE1,1)
    AFE(i,1)=AFE2(i,1)+abs(T*100-T2)*((AFE1(i,1)-AFE2(i,1))/(T2-T1));
end

% % Fit a 4th order polynomial to the interopolated discrete SA-AFE data
% Hazard_Fit_COEFF=polyfit(log(SArange),AFE,4);
% 
% % Get AFE data at all IM points
% for i=1:size(SArange,1)
%     AFE(i,1)=Hazard_Fit_COEFF(1,1)*log(SArange(i,1))^4+Hazard_Fit_COEFF(1,2)*log(SArange(i,1))^3+Hazard_Fit_COEFF(1,3)*log(SArange(i,1))^2+Hazard_Fit_COEFF(1,4)*log(SArange(i,1))^1+Hazard_Fit_COEFF(1,5);
% end

HazardCurveData(:,1)=SArange;
HazardCurveData(:,2)=exp(AFE);


