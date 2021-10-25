function Y=covid_lsq(v,t,flags,Y0,v0)

run_params;
if nargin>3
    y0 = Y0;
end

if ~isfield(flags,'cases')
    flags.cases=0;
end

switch flags.cases
    case 0
  
        params.Kc = v(1);
        params.Mc = v(2);
        params.rho0 = v(3);

        params.rhoI = 4*params.rho0;
        params.K0 = 4*params.Kc;
        params.M0 = 2*params.Mc;
        
        params.p = v0(4);
        params.rhoV0 = params.rho0/10;
        params.rhoVI = 4*params.rhoV0;
    case 1
        params.Kc = v0(1);
        params.Mc = v0(2);
        params.rho0 = v(3);

        params.rhoI = 4*params.rho0;
        params.K0 = 4*params.Kc;
        params.M0 = 2*params.Mc;
        
        params.p = v(1);
        params.rhoV0 = params.rho0/10;
        params.rhoVI = 4*params.rhoV0;
    case 2
        params.Kc = v(1);
        params.Mc = v(2);
        params.rho0 = v(3);

        params.rhoI = 4*params.rho0;
        params.K0 = 4*params.Kc;
        params.M0 = 2*params.Mc;
        
        params.p = v(4);
        params.rhoV0 = params.rho0/10;
        params.rhoVI = 4*params.rhoV0;
end

% y0 = [params.N/params.N0-v(4);0;0;0;0;0;0;0;0;0;v(4);0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];


params.beta = params.R0*params.phi*params.gamma/(params.alpha*(params.gamma+params.phi)+(1-params.alpha)*params.q*params.phi);

odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

lastwarn('');
[T,y] = ode23t(@ODEf,t,y0,odeopts,params,flags);

if ~isempty(lastwarn)
    T = t;
    y = zeros(length(t),length(y0));
end

PM = y(:,10);
IsM = y(:,14);
IaM = y(:,18);
VPM = y(:,38);
VIsM = y(:,42);
VIaM = y(:,46);
WPM = y(:,64);
WIsM = y(:,68);
WIaM = y(:,72);
M = y(:,27);

AM = PM + IsM + IaM + VPM + VIsM + VIaM + WPM + WIsM +WIaM;
VT = sum(y(:,29:80),2);
tests = params.rho0*(sum(y(:,1:10),2)+sum(y(:,15:26),2))+params.rhoI*(sum(y(:,11:14),2))+params.rhoV0*(sum(y(:,29:38),2)+sum(y(:,43:64),2)+sum(y(:,69:80),2))+params.rhoVI*(sum(y(:,39:42),2)+sum(y(:,65:68),2));

fit_dat = AM;
if flags.fit_tot
    fit_dat = [fit_dat;M];
end
if flags.fit_test
    fit_dat = [fit_dat;tests];
end

switch flags.cases
    case 0
        Y = fit_dat;
    case 1
        Y = VT;
    case 2
        Y = [fit_dat;VT];
end
end