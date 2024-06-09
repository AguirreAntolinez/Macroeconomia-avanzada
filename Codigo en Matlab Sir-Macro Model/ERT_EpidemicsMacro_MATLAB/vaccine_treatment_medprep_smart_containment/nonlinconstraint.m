function [c,ceq] = nonlinconstraint(smart_guess,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,HH,crss,nrss,Urss,phii,deltav,deltac,kappa,Uiss)

%get value of nonlinear constraint (i.e. resource constraint) must hold
%with equality (ceq)

[~,c,ceq]=getU(smart_guess,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,HH,crss,nrss,Urss,phii,deltav,deltac,kappa,Uiss);

end