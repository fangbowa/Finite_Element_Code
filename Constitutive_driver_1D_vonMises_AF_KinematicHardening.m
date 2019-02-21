%% 
% constitutive simulation using von-mises model with Armstrong-Frederick or Prager
% kinematic hardening.

clc;clear;
inc_strain=[ 0.00001*ones(1,500) -0.00001*ones(1,1000) 0.00001*ones(1,1000)];
E=1e8; % unit:Pa
ha=2e7;
Cr=1e3;
strength=1e5;

stress=zeros(1,length(inc_strain));
strain=zeros(1,length(inc_strain));
strain_pl=zeros(1,length(inc_strain));
alpha=zeros(1,length(inc_strain));
yield_tag=zeros(1,length(inc_strain));
E_ep=zeros(1,length(inc_strain));

for i=2:length(inc_strain)
    
    [stress(i), strain(i), strain_pl(i), alpha(i), yield_tag(i), E_ep(i)]= plastic_loading(stress(i-1), strain(i-1), strain_pl(i-1), alpha(i-1), inc_strain(i), yield_tag(i-1), E, ha, Cr, strength);
    
end



%%

subplot(3,1,1)
plot(strain, stress, '-*'); grid on; hold on;
xlabel('Total Strain'); ylabel('Stress (Pa)')
title({'1D von-Mises hysteresis curve (AF hardening) with',  'E=1e8 Pa, strength=1e5, ha=2e7, Cr=1e3'})

subplot(3,1,2)
plot(E_ep);
xlabel('Time step'); ylabel('Modulus (Pa)');

subplot(3,1,3)
plot(strain_pl, alpha, '-*'); grid on;
xlabel('Plastic Strain'); ylabel('Back stress (Pa)');

%%
function [stress, strain, strain_pl, alpha, yield_tag, E_ep]= plastic_loading(stress, strain, strain_pl, alpha, inc_strain, yield_tag, E, ha, Cr, strength)
    
    if yield_tag==0
        
        trialstress=stress+E*inc_strain;
        if abs(trialstress)< abs(strength+ sign(trialstress)*alpha)
            E_ep=E;
            stress=trialstress;
            strain_pl=strain_pl;
            yield_tag=0;
            alpha=alpha;
        else
            E_ep=E;
            stress=sign(trialstress)*(strength+ sign(trialstress)*alpha);
            strain_pl=strain_pl; % + (trialstress-strength+alpha(i-1))/E;
            yield_tag=1;
            alpha=alpha;
        end
    else
        L=sign(stress-alpha)*E*inc_strain;
        if L<0
           E_ep=E;
           stress=stress+E_ep*inc_strain; 
           strain_pl=strain_pl+inc_strain;
           yield_tag=0;
           alpha=alpha;
        else
           E_ep=E-E*E/(E+2/3*ha-Cr*alpha*sign(stress-alpha)); % for AF hardening
%            E_ep=E-E*E/(E+2/3*ha); % for Prager hardening
           stress=stress+E_ep*inc_strain; 
           inc_strain_pl=(1-E_ep/E)*inc_strain;
           strain_pl=strain_pl + inc_strain_pl;
           yield_tag=1;
           alpha=alpha+2/3*ha*inc_strain_pl - Cr*alpha*abs(inc_strain_pl); % for AF hardening
%            alpha=alpha+2/3*ha*inc_strain_pl; % for Prager hardening
           
        end

    end
    
    strain=strain+inc_strain;

end
