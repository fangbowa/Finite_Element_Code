%%


clc;clear;
inc_strain=[ 0.00001*ones(1,500)...
            -0.00001*ones(1,1000) 0.00001*ones(1,1000) ...
            0.00001*ones(1,500) -0.00001*ones(1,2000) 0.00001*ones(1,2000)];
E=1e8; % unit:Pa
yield_tag=0;
strength=1e5;
ha=2e7;
Cr=1e3;

stress=zeros(1,length(inc_strain));
total_strain_pl=zeros(1,length(inc_strain));
alpha=zeros(1,length(inc_strain));

for i=2:length(inc_strain)
    
    if yield_tag==0
        trialstress=stress(i-1)+E*inc_strain(i);
        if abs(trialstress)< abs(strength+ sign(trialstress)*alpha(i-1))
            E_ep=E;
            stress(i)=trialstress;
            total_strain_pl(i)=total_strain_pl(i-1);
            yield_tag=0;
            alpha(i)=alpha(i-1);
        else
            E_ep=E;
            stress(i)=sign(trialstress)*(strength+ sign(trialstress)*alpha(i-1));
            total_strain_pl(i)=total_strain_pl(i-1); % + (trialstress-strength+alpha(i-1))/E;
            yield_tag=1;
            alpha(i)=alpha(i-1);
        end
    else
        L=sign(stress(i-1)-alpha(i-1))*E*inc_strain(i);
        if L<0
           E_ep=E;
           stress(i)=stress(i-1)+E_ep*inc_strain(i); 
           total_strain_pl(i)=total_strain_pl(i-1)+inc_strain(i);
           yield_tag=0;
           alpha(i)=alpha(i-1);
        else
           E_ep=E-E*E/(E+2/3*ha-Cr*alpha(i-1)*sign(stress(i-1)-alpha(i-1)));
           stress(i)=stress(i-1)+E_ep*inc_strain(i); 
           inc_strain_pl=(1-E_ep/E)*inc_strain(i);
           total_strain_pl(i)=total_strain_pl(i-1) + inc_strain_pl;
           yield_tag=1;
           alpha(i)=alpha(i-1)+2/3*ha*inc_strain_pl - Cr*alpha(i-1)*abs(inc_strain_pl) ;
           i
        end
        
        
    end

   total_E_ep(i)=E_ep;

end



%%

subplot(3,1,1)
total_strain=cumsum(inc_strain);
plot(total_strain(1:2500), stress(1:2500), '-*'); grid on; hold on;
plot(total_strain(2500:end), stress(2500:end), '-'); grid on;
xlabel('Total Strain'); ylabel('Stress (Pa)')
legend('loop with small strain', 'loop with large strain')
title({'1D von-Mises hysteresis curve (AF hardening) with',  'E=1e8 Pa, strength=1e5, ha=2e7, Cr=1e3'})

subplot(3,1,2)

plot(total_E_ep);
xlabel('Time step'); ylabel('Modulus (Pa)');

subplot(3,1,3)

plot(total_strain_pl(1:2500), alpha(1:2500), '-*'); grid on; hold on;
plot(total_strain_pl(2500:end), alpha(2500:end)); grid on; hold on;

xlabel('Plastic Strain'); ylabel('Back stress (Pa)');


