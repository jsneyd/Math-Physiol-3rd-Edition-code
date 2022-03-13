
% This solves a closed-cell version of the calcium/secretion model, in a
% Class I version.

%% First find the steady state and plot the approach. Just for checking, basically

function main
 close all; clear all; clc; format longg; 

% The data file can be found at 
% https://github.com/jsneyd/Math-Physiol-3rd-Edition-code/blob/main/chap18_gastrointestinal/saliva_secretion_parameters.mat
 
        load('saliva_secretion_parameters.mat')            % loads all the model parameters.
        par.Ul = 10;                                      
        vplc = 0.00;        % No agonist stimulation, to calculate the steady state for resting calcium       
        M = eye(16);
        
%%%% Use these two lines if you're running this code with Octave, which doesn't deal with
## DAEs very well at all. You'll see that electroneutrality is not kept very
## well, as the cell capacitance is so high.
        M(13,13)=0.01; M(14,14)=0.01;   % 2 ODEs for Va and Vb. Capacitance is much too high for accuracy.
        options = odeset('Mass',M);

%%%% Use these two lines if you're running this code with Matlab
        %M(13,13)=0; M(14,14)=0;   % 2 DAEs for Va and Vb        
        %options = odeset('Mass',M,'RelTol', 1e-11, 'AbsTol', 1e-11);

        Nal_0     = 136.9;                                  % Na in the lumen
        Kl_0      = 6.868;                                  % K in the lumen
        HCOl_0    = 28.44;                                  % HCO in the lumen
        Hl_0      = 7.7e-5;                                 % Hl in the lumen   
        Cll_0     = Nal_0+Kl_0-HCOl_0;                      % Cl in the lumen
        w_0       = volume;                                 % cell volume
        Na_0      = 28.13;                                  % Na in the cell
        K_0       = 115.9;                                  % K in the cell
        Cl_0      = 49.76;                                  % Cl in the cell
        HCO_0     = 10.16;                                  % bicarbonate in the cell
        H_0       = 0.0001461;                              % H ions in the cell (determines the pH)
        c_init    = 0.0708;                                 % calcium concentration in the cytoplasm
        IP_0      = 0;                                      % IP3 concentration 
        HH_0      = 0.619;                                  % The IPR inactivation variable
        Va_0      = -50.71;                                 % apical membrane potential
        Vb_0      = -51.35;                                 % basal membrane potential

        IC = [Nal_0,Kl_0,Cll_0,w_0,Na_0,K_0,Cl_0,HCO_0,H_0,c_init,IP_0,HH_0,Va_0,Vb_0,HCOl_0,Hl_0];   % The initial condition

        f_secretion = @(t,x) Secretion_rhs(t,x,par,vplc,apicalarea,basalarea);
        [t,U] = ode15s(f_secretion, [0,200], IC, options);
        doplots(t,U,3,par)

%%  Now solve for a given agonist stimulation
        IC = U(end,:)';         % Start at the end point of the previous run, which will be very close to the resting steady state
        vplc = 0.0016;          % Turn on agonist stimulation for this second set of plot
        f_secretion = @(t,x) Secretion_rhs(t,x,par,vplc,apicalarea,basalarea);
        [t,U] = ode15s(f_secretion, [0,150], IC);

        doplots(t,U,4,par)
        
end


%--------------------------------------------------------------------
%%  doplots

        function doplots(t,U,figno,par)
        Nal 	= U(:,1);
        Kl 		= U(:,2);
        Cll     = U(:,3);
        w       = U(:,4);
        Na 		= U(:,5);
        K 		= U(:,6);
        Cl      = U(:,7);
        HCO     = U(:,8);
        H 		= U(:,9);
        Ca      = U(:,10);
        Ip      = U(:,11);
        Va      = U(:,13);
        Vb      = U(:,14);
        HCOl    = U(:,15);
        Hl    = U(:,16);
              
        Qa =  par.La*0.9 * ( 2 * ( Nal + Kl - Na - K - H ) - par.CO20 + par.Ul);  
        Qt =  par.Lt * ( 2 * ( Nal + Kl ) + par.Ul - par.Ie );
        Qtot=(Qa+Qt);
        
        figure
        subplot(4,5,1)
        plot(t,Qtot,'LineWidth',2)
        ylabel('fluid flow')
        mean_fluid_flow = mean(Qtot(floor(end/2):end))

        subplot(4,5,2)
        plot(t,Ca,'LineWidth',2)
        ylabel('Calcium')

        subplot(4,5,3)
        plot(t,w,'LineWidth',2)
        ylabel('cell volume')

        subplot(4,5,4)
        plot(t,Nal + Kl - Cll - HCOl,'LineWidth',2)
        ylabel('electroneutrality')
        %ylim([0 140])
        
        subplot(4,5,6)
        plot(t,Nal,'LineWidth',2)
        ylabel('Na (lumen)')
        %ylim([0 140])

        subplot(4,5,7)
        plot(t,Kl,'LineWidth',2)
        ylabel('K (lumen)')
        %ylim([0 140])

        subplot(4,5,8)
        plot(t,Cll,'LineWidth',2)
        ylabel('Cl (lumen)')
        %ylim([0 140])
        
        subplot(4,5,9)
        plot(t,HCOl,'LineWidth',2)
        ylabel('HCO (lumen)')
        
        subplot(4,5,10)
        plot(t,Nal+Kl+Cll+HCOl+par.Ul,'LineWidth',2)
        ylabel('Osmolarity (lumen)')
        %ylim([0 140])
        
        subplot(4,5,11)
        plot(t,log10(Hl*1e-3),'LineWidth',2)
        ylabel('pH (lumen)')
        %ylim([0 140])
        
        subplot(4,5,12)
        plot(t,Va,'LineWidth',2)
        ylabel('Va')

        subplot(4,5,13)
        plot(t,Vb,'LineWidth',2)
        %ylim([-60 -20])
        ylabel('Vb')
        
        subplot(4,5,14)
        plot(t,Va-Vb,'LineWidth',2)
        ylabel('Vt')
        
        subplot(4,5,15)
        plot(t,log10(H*1e-3),'LineWidth',2)
        ylabel('pH')
        %ylim([0 140])

        subplot(4,5,16)
        plot(t,K,'LineWidth',2)
        ylabel('K')
        
        subplot(4,5,17)
        plot(t,Na,'LineWidth',2)
        ylabel('Na')
        
        subplot(4,5,18)
        plot(t,Cl,'LineWidth',2)
        ylabel('Cl')
        
        subplot(4,5,19)
        plot(t,HCO,'LineWidth',2)
        ylabel('HCO')
        %ylim([-60 -20])
        
        subplot(4,5,20)
        plot(t,2*Na+2*K,'LineWidth',2)
        ylabel('Osmolarity Cell ')
        %ylim([0 140])
        
        end

%--------------------------------------------------------------------

function dx = Secretion_rhs(~,x,par,vplc,apicalarea,basalarea)

% This function contains the description of all the currents and fluxes.

Nal = x(1);
Kl = x(2);
Cll = x(3);
w = x(4);
Na = x(5);
K = x(6);
Cl = x(7);
HCO3 = x(8);
H = x(9);
Ca = x(10);
Ip = x(11);
HH = x(12);
Va = x(13);
Vb = x(14);
HCOl = x(15);
Hl = x(16);

%%%%%%%%%%%%
% IPR
ce = (1/par.gamma) * ( par.ct - Ca );
phi_c = Ca^4 / ( Ca^4 + par.K_c^4 );
phi_p = Ip^2 / ( par.K_p^2 + Ip^2);
phi_p_down = par.K_p^2 / ( par.K_p^2 + Ip^2);
H_inf = par.K_h^4 / ( par.K_h^4 + Ca^4 );
TAU = par.tau_max*par.K_tau^4/(par.K_tau^4+Ca^4);
Jserca = par.V_p*(Ca^2-par.K_bar*ce^2)/(Ca^2 + par.K_p^2);
beta = phi_p * phi_c * HH;
alpha = phi_p_down * ( 1 - phi_c * H_inf );
po = beta/(beta+par.k_beta*(beta+alpha));


%%%%%%%%%%%%%%%%%%
NaKbasalfactor = 0.7;                                                       % Fraction of NaK ATPase in the basal membrane
JNaKb = NaKbasalfactor*basalarea * par.aNaK * ( par.r * par.Ke^2 * Na^3 ...
                  / ( par.Ke^2 + par.alpha1 * Na^3 ) );
JNaKa = (1-NaKbasalfactor)*apicalarea * par.aNaK * ( par.r * Kl^2 * Na^3 ...
                  / ( Kl^2 + par.alpha1 * Na^3 ) ); 

%%%%%%%%%%%%%%%%%
VCl = par.RTF * log( Cll / Cl );    
VHCO = par.RTF * log( HCOl / HCO3 );                                      % Nernst Potentials
VKb = par.RTF * log( par.Ke / K );                                        % basal K 
VKa = par.RTF * log( Kl / K ) ;                                           % apical K 
VtNa = par.RTF * log( Nal / par.Nae );                                    % tight junction Na
VtK = par.RTF * log( Kl / par.Ke );                                       % tight junction K

Vt = Va - Vb;

%%%%%%%%%%%
% Ca2+ activated channels

totalarea = apicalarea+basalarea;
total_KCa = totalarea/1.5;                            % determined by fiddling to get nice resting volume

TMEMfiddle = 1;
KCad = 1;                                           % 1 for equal apical and basal densities
KCbd = (total_KCa - apicalarea*KCad)/basalarea;


PrCl = TMEMfiddle*1 ./ ( 1 + ( par.KCaCC ./ Ca ).^par.eta1 );  
PrKa = 1 ./ ( 1 + ( par.KCaKC ./ Ca ).^par.eta2 );    
PrKb = 1 ./ ( 1 + ( par.KCaKC ./ (Ca) ).^par.eta2 );  
PrKb = KCbd*PrKb;  
PrKa = KCad*PrKa;

JHCO = 0.3 * par.GCl * PrCl * ( Va + VHCO ) / par.F;  
JCl = par.GCl * PrCl * ( Va + VCl ) / par.F;                         % fS.micro-metres^2.mV.mol.C^-1
JKb = par.GK * PrKb * ( Vb - VKb ) / par.F;                          % basal KCa flux
JKa = 0.5*par.GK * PrKa * ( Va - VKa ) / par.F;                          % apical KCa flux

%%%%%%%%%%%%%%%%
% Tight Junction Na+ and K+ currents

JtNa = par.GtNa * par.St * ( Vt - VtNa ) / par.F;                        % fS.micro-metres^2.mV.mol.C^-1
JtK  = par.GtK * par.St * ( Vt - VtK ) / par.F;                          % fS.micro-metres^2.mV.mol.C^-1 

%%%%%%%%%%%%%%
% Water fluxes (apical, basal and tight junction)

Qa =  par.La * ( 2 * ( Nal + Kl + Hl - Na - K - H ) - par.CO20 + par.Ul );    % micro-metres^3.s^-1
Qb =  par.Lb * ( 2 * ( Na + K + H ) + par.CO20 -par.Ie);
Qt =  par.Lt * ( 2 * ( Nal + Kl + Hl ) + par.Ul - par.Ie);                    % micro-metres^3.s^-1
Qtot=(Qa+Qt);                                                            % micro-metres^3.s^-1

%%%%%%%%%%%%%
% Na+ K+ 2Cl- co-transporter (Nkcc1)

JNkcc1 = par.aNkcc1 * basalarea * ( par.a1 - par.a2 * Na * K * Cl^2 ) ...
                                             / ( par.a3 + par.a4 * Na * K * Cl^2 );     

%%%%%%%%%%%%
% (Na+)2 HCO3-/Cl- Anion exchanger (Ae4)

JAe4 = basalarea * par.G4 * ( ( par.Cle / ( par.Cle + par.KCl ) ) * ( Na / ( Na + par.KNa ) ) ...
             * ( HCO3 / ( HCO3 + par.KB ) )^2 );       

%%%%%%%%%%%%
% Na+ / H+ Anion exchanger (Nhe1)

JNhe1 = basalarea * par.G1 * ( ( par.Nae / ( par.Nae + par.KNa ) ) * ( H / ( par.KH + H ) )...
                          - ( Na / ( Na + par.KNa ) ) * ( par.He / ( par.KH + par.He ) ) ); 

%%%%%%%%%%%%%%
% Bicarbonate Buffer (Reaction Term)
% This equation is a reaction inside the cell, note how it depends on the
% cellular volume

JBB = w * par.GB * ( par.kp * par.CO20 - par.kn * HCO3 * H );
JBBA = par.wl * par.GB * ( par.kp *0.8* par.CO20 - par.kn * HCOl * Hl );

%%%%%%%%%%%%%%
% Equations

% Nal = x(1);
% Kl = x(2);
% Cll = x(3);
% w = x(4);
% Na = x(5);
% K = x(6);
% Cl = x(7);
% HCO3 = x(8);
% H = x(9);
% Ca = x(10);
% Ip = x(11);
% HH = x(12);

Jw = Qb - Qa;

dx(1) = ( JtNa - Qtot * Nal + 3*JNaKa) / par.wl;
dx(2) = ( JtK - Qtot * Kl + JKa - 2*JNaKa) / par.wl;
dx(3) = ( - JCl - Qtot * Cll ) / par.wl;
dx(4) = Jw;
dx(5) = ( JNkcc1 - 3 * (JNaKa+JNaKb) + JNhe1 - JAe4 - dx(4) * Na ) / w;
dx(6) = ( JNkcc1 + 2 * (JNaKa+JNaKb) - JKb - JKa - dx(4) * K ) / w;
dx(7) = ( 2 * JNkcc1 + JAe4 + JCl - dx(4) * Cl ) / w;
dx(8) = ( JBB - 2 * JAe4 + JHCO - dx(4) * HCO3 ) / w;
dx(9) = ( JBB - JNhe1 - dx(4) * H ) / w;
dx(10) = ((par.k_f*po )*(ce-Ca)-Jserca)-Jw*Ca / w;
dx(11) = (vplc-par.V_5K*Ip)-Jw*Ip / w;
dx(12) = (H_inf - HH) / TAU;
dx(13) = -JCl -JHCO - JNaKa - JKa - JtK - JtNa;
dx(14) = -JKb - JNaKb       + JtK + JtNa;
dx(15) = ( JBBA - JHCO - Qtot * HCOl ) / par.wl; 
dx(16) = ( JBBA - Qtot * Hl) / par.wl;

dx = dx';


end
