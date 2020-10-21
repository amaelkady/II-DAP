function [slope,force]=System_Linear (Ke, u, du)

slope=Ke;
fu=slope*du;
force=slope*u;
if force>fu
    force=0;    
end

end