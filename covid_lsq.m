function Y=covid_lsq(v,t,flags,Y0)

% persistent x;
% if isempty(x)
%     x = 1;
% else
%     x = x+1;
% end
% fprintf('Computing iteration %i...\n v=%s\n',x,num2str(v'));

% flags.model = 3;
% flags.cost_control = 0;
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

        params.rhoI = 4*v(3);

        params.K0 = 4*params.Kc;
        params.M0 = 2*params.Mc;
    case 1
        params.Kc = v(1);
        params.Mc = v(2);
        params.rho0 = v(3);

        params.rhoI = 4*v(3);

        params.K0 = 4*params.Kc;
        params.M0 = 2*params.Mc;
        
%         params.Cc = v(4);
%         params.C0 = 2*params.Cc;
%         params.numax = 0;
%         params.mumax = 100*params.mumax;
end

% y0 = [params.N/params.N0-v(4);0;0;0;0;0;0;0;0;0;v(4);0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];


params.beta = params.R0*params.phi*params.gamma/(params.alpha*(params.gamma+params.phi)+(1-params.alpha)*params.q*params.phi);

odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

[T,y] = ode23t(@ODEf,t,y0,odeopts,params,flags);

PM = y(:,10);
IsM = y(:,14);
IaM = y(:,18);
M = y(:,27);

AM = PM + IsM + IaM;

Y = [AM M];
% Y = M;
end