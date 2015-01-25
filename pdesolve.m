function [surf1,surf2,u1,u2]=pdesolve(D1, k1, D2, k2)


global roiSize;         % Size of FRAP'ed ROI (microns)
global modelLength      % Length of 1D dendrite (microns)
global deltaX           % Size of data grid cell (microns)
global simulationTime   % Simulation time (sec)
global deltaT           % Time step in the time grid
global X   % Vector [x0, x1, ..., xn] specifying the points at which a numerical solution is requested for every value in tspan ds
global T   % Vector [t0, t1,..., tf] specifying the points at which a solution is requested for every value in X
global simD2;

    
%%  Definition of the PDE system
%   pdefun() computes the terms c, f, and s (Equation 1-4). It has a form
%   [c,f,s] = pdefun(x,t,u,dudx)

    function [c,f,s] = pdefun(x,t,u,DuDx)
        c=[1; 1];               % 1 * Du_Dt 
        f=[D1; D2].*DuDx;       % Di * Du_Dx
        src=-k1*u(1)+k2*u(2);   % Source 
        s=[src;-src];
    end

%%  Initial conditions
% At t=0 there is dynamic equilibium between cytosolic and membrane HPCA
% determined by system of equations
%   u1*k1 = u2*k2 
%   u1 + u2 = 1
% This yields u1 = k2/(k1+k2), u2 = k1/(k1+k2)
    function u0 = pdeic(x)
        for i=1:length(x)
            %if abs(x(i) - modelLength/2) <= roiSize
            if abs(x(i)) <= roiSize/2
                u0(i,1)=0;
                u0(i,2)=0;
            else
                u0(i,1)=1.*k2./(k1+k2); % equilibrium conc. in membrane
                u0(i,2)=1.*k1./(k1+k2); % equilibrium conc. in cytosol
            end    
        end
    end

%%  Boundary conditions
%   Similar to initial conditions, at the ends of modeled region
%   equilibium is preserved, so 
%   ul(1) = k2/(k1+k2), ur(1) = k2/(k1+k2), - membrane 
%   ul(2) = k1/(k1+k2), ur(2) = k1/(k1+k2), - cytosol
% 
    function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
        pl=[ul(1)-1.*k2./(k1+k2);ul(2)-1.*k1./(k1+k2)];
        ql=[0;0];
        pr=[ur(1)-1.*k2./(k1+k2);ur(2)-1.*k1./(k1+k2)];
        qr=[0;0];
    end

%%  Solving PDE system and returning solution

    sol=pdepe(0,@pdefun,@pdeic,@pdebc,X,T);
    surf1=sol(:,:,1);
    surf2=sol(:,:,2);
    u1 = surf1(:,modelLength/2/deltaX);
    u2 = surf2(:,modelLength/2/deltaX);

end