function [uii,vii,aii] = ReccuranceFormulasSingleIncr(K,T,M,zeta,dt,ui,vi,agi,agii)
if dt==0
    
    uii = 0.0;
    vii = 0.0;
    aii = 0.0;   

else
    
    pi =  -M * agi;
    pii = -M * agii;

    wn=2*3.14159/T;
    wd=wn*sqrt(1-zeta^2);

    A=exp(-zeta*wn*dt)*(zeta/sqrt(1-zeta^2)*sin(wd*dt)+cos(wd*dt));
    B=exp(-zeta*wn*dt)*(sin(wd*dt)/wd);
    c1=((1-2*zeta^2)/wd/dt - zeta/sqrt(1-zeta^2))*sin(wd*dt);
    c2= (1 + 2*zeta/wn/dt)* cos(wd*dt);
    C=(2*zeta/wn/dt + exp(-zeta*wn*dt)* (c1 - c2)) /K;
    D=(1 - 2*zeta/wn/dt + exp(-zeta*wn*dt)* ((2*zeta^2-1) * sin(wd*dt)/wd/dt + 2*zeta*cos(wd*dt)/wn/dt)) /K;

    A2=-exp(-zeta*wn*dt)*(wn/sqrt(1-zeta^2)*sin(wd*dt));
    B2=exp(-zeta*wn*dt)*(cos(wd*dt) - zeta/sqrt(1-zeta^2)*sin(wd*dt));

    c1=(wn/sqrt(1-zeta^2)+ zeta/sqrt(1-zeta^2)/dt)*sin(wd*dt);
    c2=cos(wd*dt)/dt;
    C2=(-1/dt + exp(-zeta*wn*dt) * (c1+c2)) /K;
    D2=(1 - exp(-zeta*wn*dt) * (zeta/sqrt(1-zeta^2)*sin(wd*dt) + cos(wd*dt))) /K/dt;



    uii = A *ui + B *vi + C *pi + D *pii;
    vii = A2*ui + B2*vi + C2*pi + D2*pii;
    aii = (vii-vi)/dt;

end

end