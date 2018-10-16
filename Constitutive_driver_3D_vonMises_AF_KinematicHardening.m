
%%%This code is developed based on Waifah Chen's book, 
%%%Section 5.7, Section 6.5

clc;clear;


qn=1e2;  % shear strength in Pa
ha=2e7;
Cr=1000;

E_young=2e7;  % Young's modulus in Pa
poisson=0.2;
lamda=E_young*poisson/(1+poisson)/(1-2*poisson); 
miu=E_young/(1+poisson)/2;
E=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                delta_ij=0; 
                delta_kl=0;
                delta_ik=0; 
                delta_jl=0;
                delta_il=0; 
                delta_jk=0; 
                if i==j   delta_ij=1;   end
                if k==l   delta_kl=1;   end
                if i==k   delta_ik=1;   end
                if j==l   delta_jl=1;   end
                if i==l   delta_il=1;   end
                if j==k   delta_jk=1;   end
                                                                                                
                E(i,j,k,l)=lamda*delta_ij*delta_kl + miu*(delta_ik*delta_jl + delta_il*delta_jk);
            end
        end
    end
end

%pure shear strain control for the gauss point
temp_vector=[1e-5*ones(1,2000) -1e-5*ones(1,4000) 1e-5*ones(1,4000)];
for i=1:length(temp_vector)
    temp=zeros(3,3);
    stress{i}= temp;   % stress tensor
    strain{i}= temp;
    strain_pl{i}= temp;
    alpha{i}= temp;
    temp(1,2)=temp_vector(i);   %pure shear strain
    temp(2,1)=temp_vector(i);
    inc_strain{i}= temp;
end


yield_tag=0;

for i=2:length(inc_strain)
    if yield_tag==0 
        delta_stress=tensorproduct42(E, inc_strain{i});
        trialstress=stress{i-1}+ delta_stress;  
        J2=0.5*tensorproduct22(stress_to_dev(trialstress)-alpha{i-1}, stress_to_dev(trialstress)-alpha{i-1});
        f=sqrt(J2)-qn;
        if f<0
            stress{i}=trialstress; 
            strain{i}=strain{i-1}+inc_strain{i};
            strain_pl{i}=strain_pl{i-1};
            alpha{i}=alpha{i-1};
            yield_tag=0; 

        else
            
            %r = yield_function( stress{i-1}, delta_stress, alpha{i-1}, qn);
            r=0.5; % use this line in case of not converging.
            [stress{i}, strain_pl{i}, alpha{i}]=plastic_loading(E, ha, Cr, stress{i-1}+r*delta_stress, inc_strain{i}, strain_pl{i-1}, alpha{i-1});
            strain{i}=strain{i-1}+inc_strain{i};
            yield_tag=1; 
            

        end
    else
        
          J2=0.5*tensorproduct22(stress_to_dev(stress{i-1})-alpha{i-1}, stress_to_dev(stress{i-1})-alpha{i-1});
          df_dsigma=(stress_to_dev(stress{i-1})-alpha{i-1})/sqrt(J2)/2;
          L=tensorproduct22( tensorproduct42(E, inc_strain{i}), df_dsigma);
        if L<0 
            
            stress{i}=stress{i-1} + tensorproduct42(E, inc_strain{i}); 
            strain{i}=strain{i-1}+inc_strain{i};
            strain_pl{i}=strain_pl{i-1}; 
            alpha{i}=alpha{i-1};
            yield_tag=0; 

        else
          
            [stress{i}, strain_pl{i}, alpha{i}]=plastic_loading(E, ha, Cr, stress{i-1}, inc_strain{i}, strain_pl{i-1}, alpha{i-1});
            strain{i}=strain{i-1}+inc_strain{i};
            yield_tag=1; 
            
        end 
    end
    i
    
end

%%
for i=1:length(inc_strain)
    temp1=strain{i};
    temp2=stress{i};
    a(i)=temp1(1,2);
    b(i)=temp2(1,2);
end

for i=1:length(inc_strain)
    temp1=strain_pl{i};
    temp2=alpha{i};
    c(i)=temp1(1,2);
    d(i)=temp2(1,2);
end

subplot(2,1,1)
plot(a,b, '-*'); grid on;
xlabel('Shear Strain (\epsilon_{12})');
ylabel('Shear Stress (\sigma_{12})')
title({'Hystresis loop with pure shear strain control', 'E=2e7, \mu=0.2, qn=1e5 Pa, ha=2e7 Pa, Cr=1e3'})
subplot(2,1,2)
plot(c,d, '-*'); grid on;
xlabel('Plastic Shear Strain (\epsilon^{p}_{12})');
ylabel('Back Stress (\sigma_{12})')

%%
function [stress, strain_pl, alpha]=plastic_loading(E, ha, Cr, stress, inc_strain, strain_pl, alpha)

    % forward Euler method to do integration
    subincrement=100;
    for j=1:subincrement
        d_inc_strain=inc_strain/subincrement;
        J2=0.5*tensorproduct22(stress_to_dev(stress)-alpha, stress_to_dev(stress)-alpha);
        df_dsigma=(stress_to_dev(stress)-alpha)/sqrt(J2)/2;
        L=tensorproduct22( tensorproduct42(E, d_inc_strain), df_dsigma);
        h_ij = 2/3*ha*df_dsigma - Cr*alpha*sqrt(2/3);
        kappa=tensorproduct22(h_ij, df_dsigma);
        h=tensorproduct22( tensorproduct42(E, df_dsigma), df_dsigma) + kappa;
        d_lambda=L/h;
        d_inc_strain_pl=d_lambda*df_dsigma;
        d_stress= tensorproduct42(E, d_inc_strain-d_inc_strain_pl);
        stress=stress+d_stress;
        strain_pl=strain_pl+d_inc_strain_pl;
        alpha=alpha + (2/3*ha*df_dsigma - Cr*alpha*sqrt(2/3))*d_lambda;
    end


end


function C2=tensorproduct42(A4, B2)

    C2=zeros(3,3);
    for i=1:3
        for j=1:3
            temp=0;
            for k=1:3
                for l=1:3
                    temp=temp+A4(i,j,k,l)*B2(k,l);
                end
            end
            C2(i,j)=temp;
            
        end
    end
    
end

function C=tensorproduct22(A2, B2)

    C=0;
    for i=1:3
        for j=1:3           
            C=C+A2(i,j)*B2(i,j); 
        end
    end
    
end

function C2=stress_to_dev(A2)

    A2_kk=A2(1,1)+A2(2,2)+A2(3,3);
    for i=1:3
        for j=1:3     
            if i==j
               C2(i,j)=A2(i,j)-1/3*A2_kk;
            else
               C2(i,j)=A2(i,j); 
            end
        end
    end
    
end

% bisection method to find root
function r = yield_function( stress, delta_stress, alpha, qn)

    r_bound_lower=0;
    r_bound_upper=1;

    f=1;
    tol=1e-3;
    while(abs(f)>=tol)
        r=(r_bound_lower+r_bound_upper)/2;
        total_stress=stress + r*delta_stress;
        stress_dev=stress_to_dev(total_stress);

        f=sqrt( 0.5*tensorproduct22(stress_dev-alpha, stress_dev-alpha))-qn;

        if f<-tol
            r_bound_lower=r;
        elseif f>tol
            r_bound_upper=r;
        end

    end
    

end
