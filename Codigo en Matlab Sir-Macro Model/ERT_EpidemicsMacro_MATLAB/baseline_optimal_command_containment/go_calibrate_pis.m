%calibrate pis1, pis2 and pis3 using SIR model
scale1=1000000;scale2=1000;%scale pis for numerical solver
[sol,fval,exitflag]=fsolve(@calibrate_pis,[0.2;0.2;0.2],opts_fsolve,HH,i_ini,pop_ini,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,phii,crss,nrss,scale1,scale2);

if exitflag~=1
    error('Fsolve could not calibrate the SIR model');
else
    [err,pis1,pis2,pis3,RnotSIR,I,S,D,R,T] =calibrate_pis(sol,HH,i_ini,pop_ini,pir,pid,pis1_shr_target,pis2_shr_target,RplusD_target,phii,crss,nrss,scale1,scale2);
    
    disp(['Max. abs. error in calibration targets:',num2str(max(abs(err)))]);
    disp([' ']);
    pis1=sol(1)/scale1
    pis2=sol(2)/scale2
    pis3=sol(3)
    RnotSIR
end

 

% ia=2;
% ib=3;
% figure;
% subplot(ia,ib,1);
% plot(I);axis tight;
% title('I');
% subplot(ia,ib,2);
% plot(S);axis tight;
% title('S');
% subplot(ia,ib,3);
% plot(D);axis tight;
% title('D');
% subplot(ia,ib,4);
% plot(R);axis tight;
% title('R');
% subplot(ia,ib,5);
% plot(T);axis tight;
% title('T');
% suptitle('SIR Model, Calibration')
% 
% orient landscape
% print -dpdf -fillpage SIR_calibration_pis