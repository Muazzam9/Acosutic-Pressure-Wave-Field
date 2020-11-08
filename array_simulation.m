% This script solves for the normalized pressure wave field of a 2-D 
% array of rectangular elements radiating waves in a fluid using the
% MATLAB function ps_3Dv. Both time delay and apodization laws can 
% be specified for the array to steer it and focus it.


    %  ------------- give input parameters -------------------------
    
    L1 =8; % number of elements in x-direction
    %pv= 0.0024; % Lowest Voltage  10 Vp-p
    pv = 0.11;  % Highest Voltage 50Vp-p 
    
    lx = 10;   % element length in x-direction (mm)
    ly = lx;   % element length in y-direction (mm)
    f= 0.04;   % frequency (MHz)
    c = 1500;  % wave speed (m/sec)
    
    gx= ((c/(f*10^6))/2)*1000 ;      % gap length in x-direction quarter wavelength (mm)
    gy = gx; % gap length in y-direction (mm)
         
    L2 = L1;      % number of elements in y-direction
    theta =0;   % steering angle in theta direction (deg)
    phi =0;     %steering angle in phi direction (deg)
    array_length = (L1*lx)+((L1-1)*gx);
    lambda = (c/(f*10^6))*1000;
    Fl = array_length*(5/(2*lambda));  % focal distance (mm), ensure FW is 5mm
    %Fl = inf;
    disp("Focal Length (mm) = "+Fl);
    
    %weighting choices are 'rect','cos', 'Han', 'Ham', 'Blk', 'tri' 
    ampx_type ='rect';   % weighting coeffcients in x-direction
    ampy_type ='rect';   % weighting coefficients in y-direction

    % field points (x,y,z)to evaluate
    xs= linspace(-100,100, 200);
    zs= linspace(1, 300, 200);
    %ys = linspace(-100,100, 200);
    y=0;
    %z=Fl;
    %[x,y]=meshgrid(xs,ys);
    [x,z]=meshgrid(xs,zs);

    % ---------------- end input parameters ----------------------

    % calculate array pitches
    sx = lx+gx;
    sy = ly+gy;
    
    % compute centroid locations for the elements
    Nx = 1:L1;
    Ny = 1:L2;
    ex =(2*Nx -1-L1)*(sx/2);
    ey =(2*Ny -1 -L2)*(sy/2);

    % generate time delays, put in exponential 
    % and calculate amplitude weights
    td =delay_laws3D(L1,L2,sx,sy,theta,phi,Fl,c);
    delay = exp(1i.*2.*pi.*f.*td);

    Cx = discrete_windows(L1,ampx_type);
    Cy = discrete_windows(L2,ampy_type);


    % calculate normalized pressure
    p=0;
    for nn=1:L1
        for ll=1:L2
            p = p + Cx(nn)*Cy(ll)*delay(nn,ll)...
                *ps_3Dv(lx,ly,f,c,ex(nn),ey(ll),x,y,z);
                         
        end
    end
    % ---------------- outputs --------------------------
    %plot results based on specification of (x,y,z) points
    imagesc(xs,zs,abs(p)); title("2D section view");  ylabel("z in mm"); xlabel("x in mm")
    
    
    
    colormap hot
    colorbar
    
    xs= linspace(-100,100, 200);
    %zs= linspace(1, 300, 200);
    ys = linspace(-100,100, 200);
    [x,y]=meshgrid(xs,ys);
    p=0;
    for nn=1:L1
        for ll=1:L2
            p = p + Cx(nn)*Cy(ll)*delay(nn,ll)...
                *ps_3Dv(lx,ly,f,c,ex(nn),ey(ll),x,y,Fl);
                 
        end
    end
    NP = p; %Normalized pressure at center of focal point
    density = 1000;
    AP = NP.*(density*c*pv); %absolute pressure at focal point center 
    Z_0 = 1.4955*10^6 ;
    
    I = ((abs(AP).^2)./((2*Z_0)))*0.0001;
    
    figure;mesh(x,y,I);title("Intensity Profile");  ylabel("y in mm"); xlabel("x in mm");zlabel("Intensity W/cm^2")
    
    colorbar
    
 
    
  
    p=0;
    for nn=1:L1
        for ll=1:L2
            p = p + Cx(nn)*Cy(ll)*delay(nn,ll)...
                *ps_3Dv(lx,ly,f,c,ex(nn),ey(ll),0,0,Fl);
     
        end
    end
    
    
    M = max(p, [], 'all');
    phase = angle(delay);
    disp(phase)
    %disp(p)
    %disp(M);
    %disp(abs(M))
    disp("Element spacing = "+gx);
    disp("Array length = "+ array_length);

    NP = p; %Normalized pressure at center of focal point
    disp("Normalized pressure (c) = "+ NP)
    
    %Psurface = 10; %Surface pressure 1Pa
    
   
    %pv = 2.108;
    density = 1000;
    AP = NP*density*c*pv; %absolute pressure at focal point center
    disp("Absolute pressure center (Pa) = "+abs(AP));
    
    Z_0 = 1.4955*10^6 ;
    I = ((abs(AP)^2)/(2*Z_0))*0.0001; %Intensity w/cm^2 Center
   
    

    disp("Intensity center (w/cm^2) = " + I);
    

    focalpoint_width = 2*(c/40000)*1000*(Fl/array_length);
    disp("Focal point width (mm) = "+focalpoint_width);
    pe=0;
    for nn=1:L1
        for ll=1:L2
            pe = pe + Cx(nn)*Cy(ll)*delay(nn,ll)...
                *ps_3Dv(lx,ly,f,c,ex(nn),ey(ll),2.5,0,Fl);     
            

                 
        end
    end
    NPE = pe; %Normalized pressure at edge of focal point
    disp("Normalized pressure (e) = "+ NPE)
    APE = NPE*density*c*pv; %absolute pressure at focal point edge
    disp("Absolute pressure edge (Pa) = "+abs(APE));
    IE = ((abs(APE)^2)/(2*Z_0))*0.0001; %Intensity w/cm^2 Edge
    disp("Intensity edge(w/cm^2) = " + IE);
    

