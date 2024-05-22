%Script to solve TF equation using the 5 holes I.C
%from the Fig 6 in the related paper "IMEX Methods for TF equations."
%F5 will run the script without need of any other files.

dt=0.01;M1=5;iter=1;tfinal=200;
clear Energy time Mvar
M=6;
a=0;
b=M*pi;
% for proper spacing   [0 2pi]
N=256/4;
%uniform mesh thickness
h=(b-a)/N;
%(Periodic bdy conditions)
n=N;
%xgrid formtation (a b] and eventually (a b]^2
x=[a+h:h:b];[X,Y] = meshgrid(x,x);
%wave number generation (same as in 1d)
k=[[0:N/2] [-N/2+1:-1]]./((M)/2);
%now in 2-d, here k2 means k^2  and k=(k1,k2) 
 [k1x k1y]=meshgrid(k.^1,k.^1);[kx ky]=meshgrid(k.^2,k.^2);
 k2=kx+ky;k4=k2.^2;
%Initial Condition---------------------------------
U=.6-.5*exp(-.5*((X-2*pi).^2+1*(Y-1.5*pi).^2))-.5*exp(-.5*((X-4*pi).^2+1*(Y-1.5*pi).^2))+ ...
   -.5*exp(-.5*((X-2*pi).^2+1*(Y-4*pi).^2))-.5*exp(-.5*((X-4*pi).^2+1*(Y-4*pi).^2))+...
   -.5*exp(-.5*((X-3*pi).^2+1*(Y-4.5*pi).^2));

figure(1);
mesh(X,Y,U)

%parameters
epsilon=.1;eps2=epsilon^2;
M1;Bertozzinumber=max(max(U.^3));Mvar(1)=max(max(U.^3));
%The LHS, left hand side of problem----------
lhs=1+dt*M1*k4;      %Tom W.
%Left hand sidee-------------------------------

hat_U = fft2(U); 
it = 0; j=0; nn=0;   t = 0.0;
%Energy Computations
Ue1=real(ifft2(-1i*k1x.*fft2(U)));
Ue2=real(ifft2(-1i*k1y.*fft2(U)));
energy=-(eps2./(U.^2)).*((1/2)-epsilon./(3*U))+(1/2)*0.1*U.^2+(1/2)*( Ue1.^2+Ue2.^2);
Energy(1)=h*h*sum(sum(energy));
time(1)=0;
while (t <  tfinal-dt*1*0-.0000001 )

U1 = U;  %Just Eyres scheme - no initial guess      
for i=1:iter   
      MU=U1.*U1.*U1;
      phiU = (-eps2./(U1.^4)).*(3-(4*epsilon)./U1)+0.1;
      gU=MU.*phiU;       
      %-----------RHS1------------------
      %LapU
      LapU1=ifft2(-1*k2.*fft2(U1));
      %gradient (lapU1)=<lapU1_x,lapU1_y>
      lapU1_x=ifft2(1i*k1x.*fft2(LapU1));
      lapU1_y=ifft2(1i*k1y.*fft2(LapU1)); 
      %factor-in the (M1-M(u))
      lapU1_x=(M1-MU).*lapU1_x;
      lapU1_y=(M1-MU).*lapU1_y;
      rhs1=ifft2(1i*k1x.*fft2(lapU1_x))+ifft2(1i*k1y.*fft2(lapU1_y));
      rhs1=real(rhs1); 
      %-----------RHS2------------------ 
      f1=ifft2(1i*k1x.*fft2(U1));
      f2=ifft2(1i*k1y.*fft2(U1));
      %factor-in the M(u)
      f1=f1.*gU;
      f2=f2.*gU;
      
      rhs2=ifft2(1i*k1x.*fft2(f1))+ifft2(1i*k1y.*fft2(f2));
      rhs2=real(rhs2);
      %------------------------------------------------------
      RHS=rhs1+rhs2;
      hat_rhs =hat_U + dt.*fft2(RHS);
      hat_U1 = hat_rhs./lhs;
      U1 = ifft2(hat_U1);
end

U=U1;    %update 
hat_U=hat_U1;
it = it+1;
t = t+dt;
%Energy Computations
Ue1=real(ifft2(-1i*k1x.*fft2(U)));
Ue2=real(ifft2(-1i*k1y.*fft2(U)));
energy=-(eps2./(U.^2)).*((1/2)-epsilon./(3*U))+(1/2)*0.1*U.^2+(1/2)*( Ue1.^2+Ue2.^2);
Energy(it+1)=h*h*sum(sum(energy));
time(it+1)=t;
Mvar(it+1)=max(max(U.^3));
% % figure(3);                  %if monitoring the solution is needed.
% % mesh(X,Y,U)
if Energy(it+1)>Energy(it)
    Stability=0;
    break                   % This break will stop the code if energy increases - unstable
end
end  %main loop

figure(30);
mesh(X,Y,U)
title(['t_{F}=' num2str(t)] )
  U1=U;
%Energy Plotting and Reporting Stability 
figure(5)
plot(time,Energy,'.-')
axis([time(1) time(end) min(Energy) max(Energy)])
xlabel('time')
ylabel('Free Energy')
figure(6)
plot(time,Mvar,'.-')
axis([time(1) time(end) min(Mvar) max(Mvar)])
h = legend('$\bar{M}(t)$');
set(h,'interpreter','Latex','FontSize',12)
xlabel('time')
h=ylabel('$\bar{M}(t)$');
set(h,'interpreter','Latex','FontSize',12)