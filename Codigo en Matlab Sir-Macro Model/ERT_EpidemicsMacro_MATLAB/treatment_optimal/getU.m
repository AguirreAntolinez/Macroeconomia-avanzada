function [U0,err,I,S,R,D,T,Pop,cs,ns,Us,Rnot,aggC,aggH,ci,cr,ni,nr,Ui,Ur,U]=getU(muc_guess,n_vec_guess,opts_fsolve_fmincon,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uiss,HH,crss,nrss,Urss,phii,deltav,deltac,kappa)

%append guess with zeros so that muc has length HH (in case you want to
%find optimal muc with horizon less than HH)
muc=zeros(HH,1);
muc(1:numel(muc_guess))=muc_guess;
    
%create persistent variable, i.e. dont start
%solution of nonlinear equilibrium equations from scratch
%every time in the optimization of optimal containment policy
persistent guess_n n_vec
if isempty(guess_n)==1
    guess_n=n_vec_guess;
else
    guess_n=n_vec;
end;


%get allocations for ns,ni,nr given muc
try
    [n_vec,fval,exitflag]=fsolve(@get_err,guess_n,opts_fsolve_fmincon,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uiss,HH,crss,nrss,Urss,muc,phii,deltav,deltac,kappa);
                                                                
    if exitflag~=1
        error('Fsolve could not solve the model');
    end
    
    %get allocations given either exogenous or optimal path for muc at ns,ni,nr
    %solution
    [err,I,S,R,D,T,Pop,cs,ns,Us,RnotSIRmacro,aggC,aggH,ci,cr,ni,nr,Ui,Ur,U,pid_endo] = get_err(n_vec,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uiss,HH,crss,nrss,Urss,muc,phii,deltav,deltac,kappa);
                                                                                             
catch
    U=-1e10;
    guess_n=[];
end

%PV Utility at time, t=0
%Use minimizer (fmincon) to find optimal policy muc such that
%U1 is maximal -- hence flip sign
U0=-U(1);




