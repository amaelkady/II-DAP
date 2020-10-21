function [TangentK,force]=System_BiLinear (Ke, dy, Kp, fi, ui, uii)

fy=Ke*dy;
fyProj=fy-Kp*dy;

force=fi+Ke*(uii-ui);

if abs(uii) <= dy
    boundaryforceP= fyProj+Kp*uii;   
    boundaryforceN=-fyProj+Kp*uii;
end

if abs(uii) > dy
    boundaryforceP= fyProj+Kp*uii;
    boundaryforceN=-fyProj+Kp*uii;
end       

if force>=boundaryforceP
    force=boundaryforceP;
end

if force<boundaryforceN
    force=boundaryforceN;
end

TangentK=abs(force-fi)/abs(uii-ui);

if uii==ui
   TangentK=0; 
end
% 
% if TangentK>Ke
%    TangentK=Ke; 
% end
% 
% if TangentK<Kp
%    TangentK=Kp; 
% end

end