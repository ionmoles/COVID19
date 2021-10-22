function dy=ODEf(t,y,p,flags)

switch flags.model
        
    case 1
        dy = zeros(7,1);
        
        S = y(1);
        E = y(2);
        P = y(3);
        Is = y(4);
        Ia = y(5);
        Rs = y(6);
        Ra = y(7);
        
        F_S = p.alpha*p.beta*P*S+p.alpha*p.beta*Ia*S+p.beta*Is*S;
        
        dy(1) = -F_S;
        dy(2) = F_S-p.sigma*E;
        dy(3) = p.sigma*E - p.phi*P;
        dy(4) = p.q*p.phi*P - p.gamma*Is;
        dy(5) = (1-p.q)*p.phi*P - p.gamma*Ia;
        dy(6) = p.gamma*Is;
        dy(7) = p.gamma*Ia;
        
    case 2
        dy = zeros(22,1);
        
        S = y(1);
        S1 = y(2);
        S2 = y(3);
        E = y(4);
        E1 = y(5);
        E2 = y(6);
        P = y(7);
        P1 = y(8);
        P2 = y(9);
        Is = y(10);
        Is1 = y(11);
        Is2 = y(12);
        Ia = y(13);
        Ia1 = y(14);
        Ia2 = y(15);
        Rs = y(16);
        Rs1 = y(17);
        Rs2 = y(18);
        Ra = y(19);
        Ra1 = y(20);
        Ra2 = y(21);
        M = y(22);
        
        
        
        F_S = p.N0/p.N*p.beta*(p.alpha*P*S+p.alpha*Ia*S+Is*S...
            +p.alpha*p.delta*P1*S+p.alpha*p.delta*Ia1*S+p.delta*Is1*S);
        F_S1 = p.N0/p.N*p.beta*(p.alpha*p.delta*P*S1+p.alpha*p.delta*Ia*S1+p.delta*Is*S1...
            +p.alpha*p.delta^2*P1*S1+p.alpha*p.delta^2*Ia1*S1+p.delta^2*Is1*S1);
        

        mu = p.mumax*tanh(atanh(0.9)/(p.M1-p.M0)*max(M-p.M0,0));

        mu1 = 0;%0.2*mu;

        nu = -(p.numax-p.numin)*tanh(atanh(1-0.1*p.numin/(p.numax-p.numin))/(p.M1-p.M0)*max(M-p.M0,0))+p.numax;
        nu2 = nu;%1/2*nu;
        

        
        dy(1) = -F_S - mu*S + nu*S1 + (1-p.q2)*nu2*S2;
        dy(2) = -F_S1 - mu1*S1 + p.q1*mu*S - nu*S1 +p.q2*nu2*S2;
        dy(3) = (1-p.q1)*mu*S + mu1*S1 - nu2*S2;
        dy(4) = F_S - mu*E + nu*E1 + (1-p.q1)*nu2*E2 -p.sigma*E;
        dy(5) = F_S1 - mu1*E1 + p.q1*mu*E - nu*E1 + p.q1*nu2*E2 - p.sigma*E1;
        dy(6) = (1-p.q1)*mu*E + mu1*E1 - nu2*E2 - p.sigma*E2;
        dy(7) = p.sigma*E -mu*P + nu*P1 + (1-p.q1)*nu2*P2 - p.phi*P;
        dy(8) = p.sigma*E1 -mu1*P1 + p.q1*mu*P - nu*P1 + p.q1*nu2*P2 - p.phi*P1;
        dy(9) = p.sigma*E2 + (1-p.q1)*mu*P + mu1*P1 - nu2*P2 - p.phi*P2;
        dy(10) = p.q*p.phi*P - p.muI*Is - p.gamma*Is;
        dy(11) = p.q*p.phi*P1 + p.qI*p.muI*Is - p.gamma*Is1;
        dy(12) = p.q*p.phi*P2 + (1-p.qI)*p.muI*Is - p.gamma*Is2;
        dy(13) = (1-p.q)*p.phi*P - mu*Ia +nu*Ia1 + (1-p.q1)*nu2*Ia2 - p.gamma*Ia;
        dy(14) = (1-p.q)*p.phi*P1 - mu1*Ia1 + p.q1*mu*Ia - nu*Ia1 + p.q1*nu2*Ia2 - p.gamma*Ia1;
        dy(15) = (1-p.q)*p.phi*P2 + (1-p.q1)*mu*Ia + mu1*Ia1 - nu2*Ia2 - p.gamma*Ia2;
        dy(16) = p.gamma*Is;
        dy(17) = p.gamma*Is1;
        dy(18) = p.gamma*Is2;
        dy(19) = p.gamma*Ia;
        dy(20) = p.gamma*Ia1;
        dy(21) = p.gamma*Ia2;
        dy(22) = p.rho1*p.q*p.phi*(P+P1+P2) - p.muM*M;
        
    case 3
        dy = zeros(28,1);
        
        S = y(1);
        S1 = y(2);
        S2 = y(3);
        E = y(4);
        E1 = y(5);
        E2 = y(6);
        P = y(7);
        P1 = y(8);
        P2 = y(9);
        PM = y(10);
        Is = y(11);
        Is1 = y(12);
        Is2 = y(13);
        IsM = y(14);
        Ia = y(15);
        Ia1 = y(16);
        Ia2 = y(17);
        IaM = y(18);
        Rs = y(19);
        Rs1 = y(20);
        Rs2 = y(21);
        RsM = y(22);
        Ra = y(23);
        Ra1 = y(24);
        Ra2 = y(25);
        RaM = y(26);
        M = y(27);
        C = y(28);
        
        
        
%         F_S = p.N0/p.N*p.beta*(p.alpha*P*S+p.alpha*Ia*S+Is*S...
%             +p.alpha*p.delta*P1*S+p.alpha*p.delta*Ia1*S+p.delta*Is1*S);
%         F_S1 = p.N0/p.N*p.beta*(p.alpha*p.delta*P*S1+p.alpha*p.delta*Ia*S1+p.delta*Is*S1...
%             +p.alpha*p.delta^2*P1*S1+p.alpha*p.delta^2*Ia1*S1+p.delta^2*Is1*S1);

        F_S = p.N0/p.N*p.beta*S*(p.alpha*P+p.alpha*Ia+Is...
            +p.alpha*p.delta*P1+p.alpha*p.delta*Ia1+p.delta*Is1);
        F_S1 = p.N0/p.N*p.beta*S1*(p.alpha*p.delta*P+p.alpha*p.delta*Ia+p.delta*Is...
            +p.alpha*p.delta^2*P1+p.alpha*p.delta^2*Ia1+p.delta^2*Is1);
        
            KM = min(1/log(2)*max((p.rho0*(Ia+Ia1+Ia2+P+P1+P2) + p.rhoI*(Is+Is1+Is2))/M,0),1000);
            AM = PM+IsM+IaM;
%             AM = p.rho0*(P+P1+P2)+p.rhoI*(Is+Is1+Is2)+p.rho0*(Ia+Ia1+Ia2);

            mu = p.mumax*max(KM-p.Kc,0)/(max(KM-p.Kc,0)+p.K0-p.Kc)*max(AM-p.Mc,0)/(max(AM-p.Mc,0)+p.M0-p.Mc);
            mu1 = mu/2;
            
            if flags.cost_control
                nu2 = p.numax*max(C-p.Cc,0)/(max(C-p.Cc,0)+p.C0-p.Cc)*max(p.eta*p.Mc-AM,0);%/p.eta/p.Mc;
            else
                nu2 = p.numax*max(C-p.Cc,0)/(max(C-p.Cc,0)+p.C0-p.Cc);
            end
            nu = nu2/2;
            
        

        
        dy(1) = -F_S - mu*S + nu*S1 + (1-p.q2)*nu2*S2;
        dy(2) = -F_S1 - mu1*S1 + p.q1*mu*S - nu*S1 +p.q2*nu2*S2;
        dy(3) = (1-p.q1)*mu*S + mu1*S1 - nu2*S2;
        dy(4) = F_S - mu*E + nu*E1 + (1-p.q2)*nu2*E2 -p.sigma*E;
        dy(5) = F_S1 - mu1*E1 + p.q1*mu*E - nu*E1 + p.q2*nu2*E2 - p.sigma*E1;
        dy(6) = (1-p.q1)*mu*E + mu1*E1 - nu2*E2 - p.sigma*E2;
        dy(7) = p.sigma*E -mu*P + nu*P1 + (1-p.q2)*nu2*P2 - p.phi*P - p.rho0*P;
        dy(8) = p.sigma*E1 -mu1*P1 + p.q1*mu*P - nu*P1 + p.q2*nu2*P2 - p.phi*P1 - p.rho0*P1 ;
        dy(9) = p.sigma*E2 + (1-p.q1)*mu*P + mu1*P1 - nu2*P2 - p.phi*P2 - p.rho0*P2;
        dy(10) = p.rho0*(P+P1+P2) - p.phi*PM;
        dy(11) = p.q*p.phi*P - p.muI*Is - p.gamma*Is - p.rhoI*Is;
        dy(12) = p.q*p.phi*P1 + p.qI*p.muI*Is - p.gamma*Is1 - p.rhoI*Is1;
        dy(13) = p.q*p.phi*P2 + (1-p.qI)*p.muI*Is - p.gamma*Is2 - p.rhoI*Is2;
        dy(14) = p.rhoI*(Is+Is1+Is2) + p.q*p.phi*PM - p.gamma*IsM;
        dy(15) = (1-p.q)*p.phi*P - mu*Ia +nu*Ia1 + (1-p.q2)*nu2*Ia2 - p.gamma*Ia - p.rho0*Ia;
        dy(16) = (1-p.q)*p.phi*P1 - mu1*Ia1 + p.q1*mu*Ia - nu*Ia1 + p.q2*nu2*Ia2 - p.gamma*Ia1 - p.rho0*Ia1;
        dy(17) = (1-p.q)*p.phi*P2 + (1-p.q1)*mu*Ia + mu1*Ia1 - nu2*Ia2 - p.gamma*Ia2 - p.rho0*Ia2;
        dy(18) = p.rho0*(Ia+Ia1+Ia2) +(1-p.q)*p.phi*PM - p.gamma*IaM;
        dy(19) = p.gamma*Is;
        dy(20) = p.gamma*Is1;
        dy(21) = p.gamma*Is2;
        dy(22) = p.gamma*IsM;
        dy(23) = p.gamma*Ia;
        dy(24) = p.gamma*Ia1;
        dy(25) = p.gamma*Ia2;
        dy(26) = p.gamma*IsM;
        dy(27) = p.rho0*(Ia+Ia1+Ia2+P+P1+P2) + p.rhoI*(Is+Is1+Is2);
        dy(28) = p.n*((S2+E2)+(1-p.delta)*(S1+E1)) - p.muC*C;
    case 4
        dy = zeros(80,1);
        
        S = y(1);
        S1 = y(2);
        S2 = y(3);
        E = y(4);
        E1 = y(5);
        E2 = y(6);
        P = y(7);
        P1 = y(8);
        P2 = y(9);
        PM = y(10);
        Is = y(11);
        Is1 = y(12);
        Is2 = y(13);
        IsM = y(14);
        Ia = y(15);
        Ia1 = y(16);
        Ia2 = y(17);
        IaM = y(18);
        Rs = y(19);
        Rs1 = y(20);
        Rs2 = y(21);
        RsM = y(22);
        Ra = y(23);
        Ra1 = y(24);
        Ra2 = y(25);
        RaM = y(26);
        M = y(27);
        C = y(28);
        VS = y(29);
        VS1 = y(30);
        VS2 = y(31);
        VE = y(32);
        VE1 = y(33);
        VE2 = y(34);
        VP = y(35);
        VP1 = y(36);
        VP2 = y(37);
        VPM = y(38);
        VIs = y(39);
        VIs1 = y(40);
        VIs2 = y(41);
        VIsM = y(42);
        VIa = y(43);
        VIa1 = y(44);
        VIa2 = y(45);
        VIaM = y(46);
        VRs = y(47);
        VRs1 = y(48);
        VRs2 = y(49);
        VRsM = y(50);
        VRa = y(51);
        VRa1 = y(52);
        VRa2 = y(53);
        VRaM = y(54);
        WS = y(55);
        WS1 = y(56);
        WS2 = y(57);
        WE = y(58);
        WE1 = y(59);
        WE2 = y(60);
        WP = y(61);
        WP1 = y(62);
        WP2 = y(63);
        WPM = y(64);
        WIs = y(65);
        WIs1 = y(66);
        WIs2 = y(67);
        WIsM = y(68);
        WIa = y(69);
        WIa1 = y(70);
        WIa2 = y(71);
        WIaM = y(72);
        WRs = y(73);
        WRs1 = y(74);
        WRs2 = y(75);
        WRsM = y(76);
        WRa = y(77);
        WRa1 = y(78);
        WRa2 = y(79);
        WRaM = y(80);
        
        
        Infectious = p.alpha*(P+WP)+p.alpha*(Ia+WIa)+(Is+WIs)...
            +p.alpha*p.delta*(P1+WP1)+p.alpha*p.delta*(Ia1+WIa1)+p.delta*(Is1+WIs1)...
            +p.zeta*(p.alpha*VP+p.alpha*VIa+VIs+p.alpha*p.delta*VP1...
            +p.alpha*p.delta*VIa1+p.delta*VIs1);
        
        F_S = p.N0/p.N*p.beta*S*Infectious;
        F_S1 = p.N0/p.N*p.beta*p.delta*S1*Infectious;
        F_VS = p.N0/p.N*p.beta*p.epsilon*VS*Infectious;
        F_VS1 = p.N0/p.N*p.beta*p.epsilon*p.delta*VS1*Infectious;
        F_WS = p.N0/p.N*p.beta*WS*Infectious;
        F_WS1 = p.N0/p.N*p.beta*p.delta*WS1*Infectious;
        
            KM = min(1/log(2)*max((p.rho0*(Ia+Ia1+Ia2+P+P1+P2) + p.rhoI*(Is+Is1+Is2) + p.rhoV0*(VIa+WIa+VIa1+WIa1+VIa2+WIa2+VP+WP+VP1+WP1+VP2+WP2) + p.rhoVI*(VIs+WIs+VIs1+WIs1+VIs2+WIs2))/M,0),1000);
            AM = PM+IsM+IaM+VPM+VIsM+WIaM+WPM+WIsM+WIaM;
%             AM = p.rho0*(P+P1+P2)+p.rhoI*(Is+Is1+Is2)+p.rho0*(Ia+Ia1+Ia2);

            mu = max(KM-p.Kc,0)/(max(KM-p.Kc,0)+p.K0-p.Kc)*max(AM-p.Mc,0)/(max(AM-p.Mc,0)+p.M0-p.Mc);
            
            
            if flags.cost_control
                nu = max(C-p.Cc,0)/(max(C-p.Cc,0)+p.C0-p.Cc)*max(p.eta*p.Mc-AM,0);%/p.eta/p.Mc;
            else
                nu = max(C-p.Cc,0)/(max(C-p.Cc,0)+p.C0-p.Cc);
            end
            
        

        
        dy(1) = -F_S - p.mumax*mu*S + p.numax/2*nu*S1 + (1-p.q2)*p.numax*nu*S2 - p.p*S;
        dy(2) = -F_S1 - p.mumax/2*mu*S1 + p.q1*p.mumax*mu*S - p.numax/2*nu*S1 + p.q2*p.numax*nu*S2 - p.p*S1;
        dy(3) = (1-p.q1)*p.mumax*mu*S + p.mumax/2*mu*S1 - p.numax*nu*S2 - p.p*S2;
        dy(4) = F_S - p.mumax*mu*E + p.numax/2*nu*E1 + (1-p.q2)*p.numax*nu*E2 -p.sigma*E - p.p*E;
        dy(5) = F_S1 - p.mumax/2*mu*E1 + p.q1*p.mumax*mu*E - p.numax/2*nu*E1 + p.q2*p.numax*nu*E2 - p.sigma*E1 - p.p*E1;
        dy(6) = (1-p.q1)*p.mumax*mu*E + p.mumax/2*mu*E1 - p.numax*nu*E2 - p.sigma*E2 - p.p*E2;
        dy(7) = p.sigma*E -p.mumax*mu*P + p.numax/2*nu*P1 + (1-p.q2)*p.numax*nu*P2 - p.phi*P - p.rho0*P - p.p*P;
        dy(8) = p.sigma*E1 - p.mumax/2*mu*P1 + p.q1*p.mumax*mu*P - p.numax/2*nu*P1 + p.q2*p.numax*nu*P2 - p.phi*P1 - p.rho0*P1 - p.p*P1;
        dy(9) = p.sigma*E2 + (1-p.q1)*p.mumax*mu*P + p.mumax/2*mu*P1 - p.numax*nu*P2 - p.phi*P2 - p.rho0*P2 - p.p*P2;
        dy(10) = p.rho0*(P+P1+P2) - p.phi*PM;
        dy(11) = p.q*p.phi*P - p.muI*Is - p.gamma*Is - p.rhoI*Is;
        dy(12) = p.q*p.phi*P1 + p.qI*p.muI*Is - p.gamma*Is1 - p.rhoI*Is1;
        dy(13) = p.q*p.phi*P2 + (1-p.qI)*p.muI*Is - p.gamma*Is2 - p.rhoI*Is2;
        dy(14) = p.rhoI*(Is+Is1+Is2) + p.q*p.phi*PM - p.gamma*IsM;
        dy(15) = (1-p.q)*p.phi*P - p.mumax*mu*Ia + p.numax/2*nu*Ia1 + (1-p.q2)*p.numax*nu*Ia2 - p.gamma*Ia - p.rho0*Ia - p.p*Ia;
        dy(16) = (1-p.q)*p.phi*P1 - p.mumax/2*mu*Ia1 + p.q1*p.mumax*mu*Ia - p.numax/2*nu*Ia1 + p.q2*p.numax*nu*Ia2 - p.gamma*Ia1 - p.rho0*Ia1 - p.p*Ia1;
        dy(17) = (1-p.q)*p.phi*P2 + (1-p.q1)*p.mumax*mu*Ia + p.mumax/2*mu*Ia1 - p.numax*nu*Ia2 - p.gamma*Ia2 - p.rho0*Ia2 - p.p*Ia2;
        dy(18) = p.rho0*(Ia+Ia1+Ia2) +(1-p.q)*p.phi*PM - p.gamma*IaM;
        dy(19) = p.gamma*Is - p.p*Rs;
        dy(20) = p.gamma*Is1 - p.p*Rs1;
        dy(21) = p.gamma*Is2 - p.p*Rs2;
        dy(22) = p.gamma*IsM - p.p*RsM;
        dy(23) = p.gamma*Ia - p.p*Ra;
        dy(24) = p.gamma*Ia1 - p.p*Ra1;
        dy(25) = p.gamma*Ia2 - p.p*Ra2;
        dy(26) = p.gamma*IaM - p.p*RaM;
        dy(27) = p.rho0*(Ia+Ia1+Ia2+P+P1+P2) + p.rhoI*(Is+Is1+Is2) + p.rhoV0*(VIa+WIa+VIa1+WIa1+VIa2+WIa2+VP+WP+VP1+WP1+VP2+WP2) + p.rhoVI*(VIs+WIs+VIs1+WIs1+VIs2+WIs2);
        dy(28) = p.n*((S2+E2+VS2+VE2+WS2+WE2)+(1-p.delta)*(S1+E1+VS1+VE1+WE1+WE2)) - p.muC*C;
        
        %Vaccination compartments start here
        dy(29) = -F_VS - p.mumaxV*mu*VS + p.numaxV/2*nu*VS1 + (1-p.q2)*p.numaxV*nu*VS2 + p.qVS*p.p*S;
        dy(30) = -F_VS1 - p.mumaxV/2*mu*VS1 + p.q1*p.mumaxV*mu*VS - p.numaxV/2*nu*VS1 + p.q2*p.numaxV*nu*VS2 + p.qVS*p.p*S1;
        dy(31) = (1-p.q1)*p.mumaxV*mu*VS + p.mumaxV/2*mu*VS1 - p.numaxV*nu*VS2 + p.qVS*p.p*S2;
        dy(32) = F_VS - p.mumaxV*mu*VE + p.numaxV/2*nu*VE1 + (1-p.q2)*p.numaxV*nu*VE2 -p.sigma*VE + p.qVE*p.p*E;
        dy(33) = F_VS1 - p.mumaxV/2*mu*VE1 + p.q1*p.mumaxV*mu*VE - p.numaxV/2*nu*VE1 + p.q2*p.numaxV*nu*VE2 - p.sigma*VE1 + p.qVE*p.p*E1;
        dy(34) = (1-p.q1)*p.mumaxV*mu*VE + p.mumaxV/2*mu*VE1 - p.numaxV*nu*VE2 - p.sigma*VE2 + p.qVE*p.p*E2;
        dy(35) = p.sigma*VE -p.mumaxV*mu*VP + p.numaxV/2*nu*VP1 + (1-p.q2)*p.numaxV*nu*VP2 - p.phi*VP - p.rhoV0*VP + p.qVP*p.p*P;
        dy(36) = p.sigma*VE1 - p.mumaxV/2*mu*VP1 + p.q1*p.mumaxV*mu*VP - p.numaxV/2*nu*VP1 + p.q2*p.numaxV*nu*VP2 - p.phi*VP1 - p.rhoV0*VP1 + p.qVP*p.p*P1;
        dy(37) = p.sigma*VE2 + (1-p.q1)*p.mumaxV*mu*VP + p.mumaxV/2*mu*VP1 - p.numaxV*nu*VP2 - p.phi*VP2 - p.rhoV0*VP2 + p.qVP*p.p*P2;
        dy(38) = p.rhoV0*(VP+VP1+VP2) - p.phi*VPM;
        dy(39) = p.qv*p.phi*VP - p.muI*VIs - p.gamma*VIs - p.rhoVI*VIs;
        dy(40) = p.qv*p.phi*VP1 + p.qI*p.muI*VIs - p.gamma*VIs1 - p.rhoVI*VIs1;
        dy(41) = p.qv*p.phi*VP2 + (1-p.qI)*p.muI*VIs - p.gamma*VIs2 - p.rhoVI*VIs2;
        dy(42) = p.rhoVI*(VIs+VIs1+VIs2) + p.qv*p.phi*VPM - p.gamma*VIsM;
        dy(43) = (1-p.qv)*p.phi*VP - p.mumaxV*mu*VIa + p.numaxV/2*nu*VIa1 + (1-p.q2)*p.numaxV*nu*VIa2 - p.gamma*VIa - p.rhoV0*VIa + p.qVIa*p.p*Ia;
        dy(44) = (1-p.qv)*p.phi*VP1 - p.mumaxV/2*mu*VIa1 + p.q1*p.mumaxV*mu*VIa - p.numaxV/2*nu*VIa1 + p.q2*p.numaxV*nu*VIa2 - p.gamma*VIa1 - p.rhoV0*VIa1 + p.qVIa*p.p*Ia1;
        dy(45) = (1-p.qv)*p.phi*VP2 + (1-p.q1)*p.mumaxV*mu*VIa + p.mumaxV/2*mu*VIa1 - p.numaxV*nu*VIa2 - p.gamma*VIa2 - p.rhoV0*VIa2 + p.qVIa*p.p*Ia2;
        dy(46) = p.rhoV0*(VIa+VIa1+VIa2) +(1-p.qv)*p.phi*VPM - p.gamma*VIaM;
        dy(47) = p.gamma*VIs + p.qVR*p.p*Rs;
        dy(48) = p.gamma*VIs1 + p.qVR*p.p*Rs1;
        dy(49) = p.gamma*VIs2 + p.qVR*p.p*Rs2;
        dy(50) = p.gamma*VIsM + p.qVR*p.p*RsM;
        dy(51) = p.gamma*VIa + p.qVR*p.p*Ra;
        dy(52) = p.gamma*VIa1 + p.qVR*p.p*Ra1;
        dy(53) = p.gamma*VIa2 + p.qVR*p.p*Ra2;
        dy(54) = p.gamma*VIaM + p.qVR*p.p*RaM;     
        
        dy(55) = -F_WS - p.mumaxV*mu*WS + p.numaxV/2*nu*WS1 + (1-p.q2)*p.numaxV*nu*WS2 + (1-p.qVS)*p.p*S;
        dy(56) = -F_WS1 - p.mumaxV/2*mu*WS1 + p.q1*p.mumaxV*mu*WS - p.numaxV/2*nu*WS1 + p.q2*p.numaxV*nu*WS2 + (1-p.qVS)*p.p*S1;
        dy(57) = (1-p.q1)*p.mumaxV*mu*WS + p.mumaxV/2*mu*WS1 - p.numaxV*nu*WS2 + (1-p.qVS)*p.p*S2;
        dy(58) = F_WS - p.mumaxV*mu*WE + p.numaxV/2*nu*WE1 + (1-p.q2)*p.numaxV*nu*WE2 - p.sigma*WE + (1-p.qVE)*p.p*E;
        dy(59) = F_WS1 - p.mumaxV/2*mu*WE1 + p.q1*p.mumaxV*mu*WE - p.numaxV/2*nu*WE1 + p.q2*p.numaxV*nu*WE2 - p.sigma*WE1 + (1-p.qVE)*p.p*E1;
        dy(60) = (1-p.q1)*p.mumaxV*mu*WE + p.mumaxV/2*mu*WE1 - p.numaxV*nu*WE2 - p.sigma*WE2 + (1-p.qVE)*p.p*E2;
        dy(61) = p.sigma*WE -p.mumaxV*mu*WP + p.numaxV/2*nu*WP1 + (1-p.q2)*p.numaxV*nu*WP2 - p.phi*WP - p.rhoV0*WP + (1-p.qVP)*p.p*P;
        dy(62) = p.sigma*WE1 - p.mumaxV/2*mu*WP1 + p.q1*p.mumaxV*mu*WP - p.numaxV/2*nu*WP1 + p.q2*p.numaxV*nu*WP2 - p.phi*WP1 - p.rhoV0*WP1 + (1-p.qVP)*p.p*P1;
        dy(63) = p.sigma*WE2 + (1-p.q1)*p.mumaxV*mu*WP + p.mumaxV/2*mu*WP1 - p.numaxV*nu*WP2 - p.phi*WP2 - p.rhoV0*WP2 + (1-p.qVP)*p.p*P2;
        dy(64) = p.rhoV0*(WP+WP1+WP2) - p.phi*WPM;
        dy(65) = p.q*p.phi*WP - p.muI*WIs - p.gamma*WIs - p.rhoVI*WIs;
        dy(66) = p.q*p.phi*WP1 + p.qI*p.muI*WIs - p.gamma*WIs1 - p.rhoVI*WIs1;
        dy(67) = p.q*p.phi*WP2 + (1-p.qI)*p.muI*WIs - p.gamma*WIs2 - p.rhoVI*WIs2;
        dy(68) = p.rhoVI*(WIs+WIs1+WIs2) + p.q*p.phi*WPM - p.gamma*WIsM;
        dy(69) = (1-p.q)*p.phi*WP - p.mumaxV*mu*WIa + p.numaxV/2*nu*WIa1 + (1-p.q2)*p.numaxV*nu*WIa2 - p.gamma*WIa - p.rhoV0*WIa + (1-p.qVIa)*p.p*Ia;
        dy(70) = (1-p.q)*p.phi*WP1 - p.mumaxV/2*mu*WIa1 + p.q1*p.mumaxV*mu*WIa - p.numaxV/2*nu*WIa1 + p.q2*p.numaxV*nu*WIa2 - p.gamma*WIa1 - p.rhoV0*WIa1 + (1-p.qVIa)*p.p*Ia1;
        dy(71) = (1-p.q)*p.phi*WP2 + (1-p.q1)*p.mumaxV*mu*WIa + p.mumaxV/2*mu*WIa1 - p.numaxV*nu*WIa2 - p.gamma*WIa2 - p.rhoV0*WIa2 + (1-p.qVIa)*p.p*Ia2;
        dy(72) = p.rhoV0*(WIa+WIa1+WIa2) +(1-p.q)*p.phi*WPM - p.gamma*WIaM;
        dy(73) = p.gamma*WIs + (1-p.qVR)*p.p*Rs;
        dy(74) = p.gamma*WIs1 + (1-p.qVR)*p.p*Rs1;
        dy(75) = p.gamma*WIs2 + (1-p.qVR)*p.p*Rs2;
        dy(76) = p.gamma*WIsM + (1-p.qVR)*p.p*RsM;
        dy(77) = p.gamma*WIa + (1-p.qVR)*p.p*Ra;
        dy(78) = p.gamma*WIa1 + (1-p.qVR)*p.p*Ra1;
        dy(79) = p.gamma*WIa2 + (1-p.qVR)*p.p*Ra2;
        dy(80) = p.gamma*WIaM + (1-p.qVR)*p.p*RaM;     
        
end


