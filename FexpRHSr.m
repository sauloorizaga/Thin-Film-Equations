function [RHS] = FexpRHSr(epsilon,eps2,M1,k1x,k1y,k2,U1)

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
      %U1(U1<=epsilon/100)=epsilon/100;   
 
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
      %------------------------------------------------------------------
      
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

end
