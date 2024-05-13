%This code will serve as a tool to determine
%energy-stable solutions.

function [Stability]=BertoziMethod_Stability(dt,M1,iter,tfinal)

clear Energy time Mvar GUvar
clc;

Stability=1;

M=12;
a=0;
b=M*pi;
%number of grid points N
N=256;
%uniform mesh thickness
h=(b-a)/N;
%(Periodic bdy conditions)
n=N;

%xgrid formtation (a b] and eventually (a b]^2
x=[a+h:h:b];[X,Y] = meshgrid(x,x);

%wave number generation (same as in 1d)
k=[[0:N/2] [-N/2+1:-1]]./((M+4*0)/2);

%now in 2-d, here k2 means k^2  and k=(k1,k2) 
 [k1x k1y]=meshgrid(k.^1,k.^1);
 
 [kx ky]=meshgrid(k.^2,k.^2);
 k2=kx+ky;
 k4=k2.^2;

%Initial Condition---------------------------------
load('ICactive.mat','U');
%--------------------------------------------------

figure(1);
mesh(X,Y,U)
ax = gca; 
ax.FontSize = 14;
colormap('jet')

%parameters
epsilon=.1*1;
eps2=epsilon^2;

M1;
Bertozzinumber=max(max(U.^3));
Mvar(1)=max(max(U.^3));
GUvar(1)= max(max(-epsilon*epsilon*(3.0.*log(U)+4.0*epsilon./U)+0.25.*(0.1).*U.^4));

%The LHS, left hand side of problem----------
lhs=1+dt*M1*k4;      %Tom W.
%lhs=1+dt*M1*k4*eps2;   %A.Bertozzi
%Left hand sidee-------------------------------


hat_U = fft2(U); 
it = 0; j=0; nn=0;   t = 0.0;

%Energy Computations
Ue1=real(ifft2(-1i*k1x.*fft2(U)));
Ue2=real(ifft2(-1i*k1y.*fft2(U)));
energy=-(eps2./(U.^2)).*((1/2)-epsilon./(3*U))+(1/2)*0.1*U.^2+(1/2)*( Ue1.^2+Ue2.^2);
%%energy=(1/2)*( Ue1.^2+Ue2.^2);     %Lubrication Only
Energy(1)=h*h*sum(sum(energy));
time(1)=0;


while (t <  tfinal-dt*1*0-.0000001 )

U1 = U;  %Just Eyres scheme - no initial guess
             
for i=1:iter
      RHS=FexpRHSr(epsilon,eps2,M1,k1x,k1y,k2,U1);
      
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
%%energy=(1/2)*( Ue1.^2+Ue2.^2);     %Lubrication Only
Energy(it+1)=h*h*sum(sum(energy));
time(it+1)=t;

Mvar(it+1)=max(max(U.^3));
GUvar(it+1)= max(max(-eps*eps*(3.0.*log(U)+4.0*eps./U)+0.25.*(0.1).*U.^4));
GUvar(it+1)= max(max(-epsilon*epsilon*(3.0.*log(U)+4.0*epsilon./U)+0.25*(0.1).*U.^4));

% % figure(3);                  %if monitoring the solution is needed.
% % mesh(X,Y,U)

if Energy(it+1)>Energy(it)
    Stability=0;
    break                   % This break will stop the code if energy increases - unstable
end
   
% 
% if M1>.3 && dt<20/2^6
%     Stability=1;
%     break                   % This break code for already known stanle regions
% end
% 
% if M1>.28 && dt<20/2^12
%     Stability=1;
%     break                   % This break code for already known stanle regions
% end
% 
% if M1<.22 && dt>20/2^14
%     Stability=0;
%     break                   % This break code for already known stanle regions
% end


end  %main loop

figure(30);
mesh(X,Y,U)
ax = gca; 
ax.FontSize = 14;
colormap('jet')
%title(['t_{F}=' num2str(t)] )


  U1=U;
  load('Uexact','U');
  error=h*h.*sum(sum(abs(U1-U)));          %L1 error
 
 figure(31);
mesh(X,Y,U)
ax = gca; 
ax.FontSize = 14;
colormap('jet')
%title(['t_{F}=' num2str(t)] )

 %Energy Plotting and Reporting Stability 
 % To be commented out for efficiency
% % figure(5)
% % plot(time,Energy,'.-')
% % axis([time(1) time(end) min(Energy) max(Energy)])
% % xlabel('time')
% % ylabel('Free Energy')
% % 
% % figure(6)
% % %plot(time,Mvar,'o-',time,GUvar,'o-')
% % plot(time,Mvar,'.-')
% % axis([time(1) time(end) min(Mvar) max(Mvar)])
% % h = legend('$\bar{M}(t)$');
% % set(h,'interpreter','Latex','FontSize',12)
% % %legend('$\bar M(t)$')
% % xlabel('time')
% % %ylabel('M_1')
% % h=ylabel('$\bar{M}(t)$');
% % set(h,'interpreter','Latex','FontSize',12)
% % 
% % figure(7)
% % plot(time,GUvar,'.-')
% % axis([time(1) time(end) min(GUvar) max(GUvar)])
% % h = legend('$\bar{G}(t)$');
% % set(h,'interpreter','Latex','FontSize',12)
% % %legend('$\bar M(t)$')
% % xlabel('time')
% % %ylabel('G_1')
% % h=ylabel('$\bar{G}(t)$');
% % set(h,'interpreter','Latex','FontSize',12)
end