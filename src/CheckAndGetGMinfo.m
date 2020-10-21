% GMfilePath    Path to ground motion file
% GMfilename    Name of ground motion file
% GMx           Structured variable containing ground motion data
% GMx.Time      Ground motion time history
% GMx.GA        Ground motion acceleration history
% GMx.npoints   Ground motion number of points
% GMx.duration  Ground motion duration
% GMx.pga       Ground motion peak acceleration
% GMx.dt        Ground motion time step

function    [GMx, errormsg, inputtype]=CheckAndGetGMinfo(GMfilePath, GMfilename)

    MainDir=pwd;
	cd (GMfilePath)
    data=importdata(GMfilename);

    if size(data,2)==2
        
        inputtype=1;
        GMx.name={GMfilename};
        GMx.Time=data(:,1);
        GMx.GA=data(:,2);
        GMx.npoints=size(data,1);
        GMx.duration=data(end,1);
        GMx.pga=max(abs(data(:,2)));
        GMx.dt=data(end,1)-data(end-1,1);

        x=diff(diff(data(:,1)));
        iDX=find(x>1^-10);
        if size(iDX,1)~=0
            errormsg=1;
        else
            errormsg=0;
        end
    
    elseif size(data,2)==1
        
        errormsg=0;
        inputtype=2;
        GMx.name={GMfilename};
        GMx.GA=data(:,1);
        GMx.npoints=size(data,1);
        GMx.pga=max(abs(data(:,1)));
        
    elseif size(data,2)>2
        
        counter=0;
        for i=1:size(data,1)
           for j=1:size(data,2)
                counter=counter+1;
                GMx.GA(counter,1)=data(i,j); 
           end
        end
        errormsg=0;
        inputtype=2;
        GMx.name={GMfilename};
        GMx.npoints=counter;
        GMx.pga=max(abs(GMx.GA(:,1))); 
        
    end
    
    cd (MainDir)
end