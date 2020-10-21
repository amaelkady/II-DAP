function [uii,vii,aii, ui_1] = CentralDifferenceSingleIncrNL(C,M,dt,ui_1,ui,agii,fii)

if dt==0
    uii = 0.0;
    vii = 0.0;
    aii = 0.0;   
else
    pii = -M * agii;

    A     = M/dt^2 - C/2/dt;
    B     = 2*M/dt^2;
    k_eff = M/dt^2 + C/2/dt;

    p_eff = pii - A*ui_1 + B*ui -fii;

    uii = p_eff/k_eff;
    vii = (uii-ui_1)/2/dt;
    aii = (uii-2*ui+ui_1)/dt^2;
end

end