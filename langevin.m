function sol=langevin(kernel, steps, par, u)

% Simulation of collective abraison at master equation level
% This code uses the ChebFun Toolbox (https://www.chebfun.org/)
% Input:
%   kernel  :function handle with two arguments
%   steps   :number of timesteps
%   par     :array of control parameters, see the explanations below
%   u       :initial distribution (any continuous chebfun)
% Output:
%   sol     :result structure
% 
% Example:
% sol=langevin(@(x,y)-x-y,1000)
%
%Andras A. Sipos
%2021.01.01.

%initialization
if nargin<3 
    dt=0.01;            %Time step
    ac=1;               %Accuracy parameter
    dp=0.000;           %diffusion parameter (suggested:10^-3 or 0)
    par=[dt,ac,dp];        
else
    dt=par(1);
    ac=par(2);
    dp=par(3);
end

if isa(kernel,'double')==1  %initation of the compound kernel
    par(4)=kernel;
    kernel=@(x,z)-x^(1+kernel)*z^(1-kernel)/(x+z);    
end

%initial distribution
if nargin<4
    xmin=0.001;  %cut-off at small size
    xmax=5.0;    %cut-off at large size
    x = chebfun('x',[xmin xmax]);    
    %s2=10;
    %u=chebfun(sqrt(s2/pi)*exp(-s2*(x-1)^2),[xmin,xmax],'splitting','off')
    %u=chebfun(0.5*sqrt(s2/pi)*exp(-s2*(x-1)^2)+0.5*sqrt(s2/pi)*exp(-s2*(x-3)^2),[xmin,xmax],'splitting','off');  
    mu=-1/32;
    sig=1/4;
    u=chebfun(1/sig/sqrt(2*pi)/x.*exp(-(log(x)-mu).^2/2/sig^2),[xmin,xmax],'splitting','off') %lognormal initial 
end

%kernel initialization
ker =chebfun2(@(x,z) kernel(x,z),[xmin,xmax,xmin,xmax],'vectorize');
sol(1).u=u;
E(1)=sum(sol(1).u.*x);
V(1)=sum(sol(1).u.*x.*x)-E(1)^2;
t(1)=0;

%a linear operator
L = chebop(xmin,xmax);

%evolution
w=2;
if dp>0  %evolution with diffusion
    L.lbc=@(v) diff(v,2);
    L.rbc=@(v) diff(v,2);
    for k=1:steps
        for l=1:ac
            D1=fred(ker,u)';
            L.op=@(x,v) -diff(D1.*v,1)+dp*diff(v,2);
            u=expm(L, dt, u);
        end
            sol(w).u=u;
            E(w)=sum(sol(w).u.*x);
            V(w)=sum(sol(w).u.*x.*x)-E(w).^2;
            t(w)=k*dt;
            w=w+1      
    end
else    %evolution without diffusion
    L.rbc=0;
    for k=1:steps
        for l=1:ac
            D1=fred(ker,u)';
            L.op=@(x,v) -diff(D1.*v,1);
            u=expm(L, dt, u);
        end
            sol(w).u=u;
            E(w)=sum(sol(w).u.*x);
            V(w)=sum(sol(w).u.*x.*x)-E(w).^2;
            t(w)=k*dt;
            w=w+1      
    end
end    

%writing results
sol(1).E=E;
sol(1).V=V;
sol(1).t=t;
sol(1).par=par;
sol(1).r=par(end);