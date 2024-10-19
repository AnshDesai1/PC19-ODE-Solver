function [u,t,km,counter]=OneStep(f,u,t,tol1,H,counter)
done=0;
l=0;
maxitn=100;
while ~done
   % Compute the slopes:
   if l==0
       k1=f(t(1),u(:,1));
       counter=counter+1;
   end
   k2=f(t(1)+H,u(:,1)+H*k1);
   counter=counter+1;
   % Compute the first-order estimate:
   u1=u(:,1)+H*k1;
   % Compute the second-order estimate:
   u2=u(:,1)+H*(0.5*k1+0.5*k2);
   % Estimate the local truncation error:
   lte=norm(u2-u1);
   % Decide how to proceed:
   l=l+1;
   if lte>=tol1 && l<maxitn
      % Reject this step, reduce h, and try again:
       H=H*sqrt(tol1/(2*lte));
   elseif l==maxitn
      error('Maximum number of step-size reductions reached in initialization.')
   else
      % Accept this step:
      t(2)=t(1)+H;
      u(:,2)=u2;
      km(:,1)=k1;
      done=1;
   end
end
end