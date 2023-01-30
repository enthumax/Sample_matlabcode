function [ ncoeff,fcoeff,func, nv,vvec,PV_norm,vglobal1 ] = normlaizationwithfilter_costfn2_updated_area( X1,Y,alphayesno,s1,s2,dsprange,prefdsp,B)
%Divisive normalization model with a gaussian filter for s1 and s2
%   Detailed explanation goes here
% y = normpdf(x,mu,sigma) 

[vglobal1, area, vec]=vglobal_area(s1,s2,dsprange,prefdsp);  %or vglobal function


if alphayesno=="No"
beta0 = [1 0 0]; % beta(1) = n, beta(2) = sigma, beta(3) = C
        A = [];
        b = [];
        Aeq = []; 
        beq = [];
        %look for filter width v within disparity range -1.6 to 1.6 deg
        lb = [0 0 0 ];   % changed the sigma lower bound from -inf to 0 on 1/9/2020
        ub = [100 5000 100];    % changed these upper bounds on 1/9/2020
        lb1 = 0.01;   % changed the sigma lower bound from -inf to 0 on 1/9/2020
        ub1 = 15; %15;    % If we assume 3.2 is the max full width at half max v (or std)=3.2/2.355
        
        s=0.1;
        i=0;
        tol=0.08; %tolerance to prefdsp value incase there is no exact match for prefdsp in dsprange
        %costfun=
        for step=[0 lb1:s:ub1 100]
       i=i+1;
       if step==0
           fun = @(w)sum(((B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^w(1).*X1(:,1) +((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^w(1)).*X1(:,2))./ ...
                      (B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^w(1) + ((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^w(1)) + w(2)) + w(3) - Y).^2);
        [beta, func(i)] = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        
        nv(i)=beta(1);
        ncoeff(i)=sum(s1(ismembertol(dsprange,prefdsp,tol))).^beta(1);
        fcoeff(i)=sum(s2(ismembertol(dsprange,prefdsp,tol))).^beta(1);
        predD60_N_c_sort = (B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^beta(1).*X1(:,1) + ((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^beta(1)).*X1(:,2))./ ...
                (B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^beta(1) + ((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^beta(1)) + beta(2)) + beta(3);
            SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_norm(i) = 100*(1 - SSR/SST);
         else
           PV_norm(i) = 0;
         end;
           
       elseif step==100
           fun = @(w)sum(((B*(sum(s1)).^w(1).*X1(:,1) +((sum(s2)).^w(1)).*X1(:,2))./ ...
                      (B*(sum(s1)).^w(1) + ((sum(s2)).^w(1)) + w(2)) + w(3) - Y).^2);
        [beta, func(i)] = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        
        nv(i)=beta(1);
        ncoeff(i)=sum(s1).^beta(1);
        fcoeff(i)=sum(s2).^beta(1);
        predD60_N_c_sort = (B*(sum(s1)).^beta(1).*X1(:,1) + ((sum(s2)).^beta(1)).*X1(:,2))./ ...
                (B*(sum(s1)).^beta(1) + ((sum(s2)).^beta(1)) + beta(2)) + beta(3);
            SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_norm(i) = 100*(1 - SSR/SST);
         else
           PV_norm(i) = 0;
         end;
       else
       fun = @(w)sum(((B*(sum(s1.*normpdf(dsprange,prefdsp,step))).^w(1).*X1(:,1) +((sum(s2.*normpdf(dsprange,prefdsp,step))).^w(1)).*X1(:,2))./ ...
                      (B*(sum(s1.*normpdf(dsprange,prefdsp,step))).^w(1) + ((sum(s2.*normpdf(dsprange,prefdsp,step))).^w(1)) + w(2)) + w(3) - Y).^2);
        [beta, func(i)] = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        
        nv(i)=beta(1);
        ncoeff(i)=sum(s1.*normpdf(dsprange,prefdsp,step)).^beta(1);
        fcoeff(i)=sum(s2.*normpdf(dsprange,prefdsp,step)).^beta(1);
        
        predD60_N_c_sort = (B*(sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1).*X1(:,1) + ((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)).*X1(:,2))./ ...
                (B*(sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1) + ((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)) + beta(2)) + beta(3);
            SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_norm(i) = 100*(1 - SSR/SST);
         else
           PV_norm(i) = 0;
         end;
     
%         %add four conditions up
%         for step=lb:s:ub
%        i=i+1;
%        func(i) = sum((((sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1).*X1(:,1) +((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)).*X1(:,2))./ ...
%                       ((sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1) + ((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)) + beta(2)) + beta(3) - Y).^2);
%         
       end
        end
        vvec=[0 lb1:s:ub1 ub1+1];
        
        
        
elseif alphayesno=="Yes"
    
    %find the global value for each prefdsp i.e.,
    %(sum(s1.*normpdf(dsprange,prefdsp,v))-sum(s1)).^2 is minimized
    
    A = [];
        b = [];
        Aeq = []; 
        beq = [];
    
    beta0 = [1 0 0 0.5]; % beta(1) = n, beta(2) = sigma, beta(3) = C
        
        %look for filter width v within disparity range -1.6 to 1.6 deg
        lb = [0 0 0 0.01];   %alpha 0.1 to10
        ub = [100 5000 100 100];    % changed these upper bounds on 1/9/2020
        lb1 = min(area);   % changed the sigma lower bound from -inf to 0 on 1/9/2020
        ub1 = max(area); %15;    % If we assume 3.2 is the max full width at half max v (or std)=3.2/2.355
        
        s=(max(area)-min(area))/20;
        i=0;
        tol=0.08; %tolerance to prefdsp value incase there is no exact match for prefdsp in dsprange
        %costfun=
        for step=[lb1:s:ub1]
       i=i+1;
         
         tol1=0.01;
         v=min(vec(ismembertol(area,step,tol1)));
       if step==0
           fun = @(w)sum(((B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^w(1).*X1(:,1) +((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^w(1)).*X1(:,2))./ ...
                      (B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^w(1) + w(4).*((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^w(1)) + w(2)) + w(3) - Y).^2);
        [beta, func(i)] = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        
        nv(i)=beta(1);
        ncoeff(i)=sum(s1(ismembertol(dsprange,prefdsp,tol))).^beta(1);
        fcoeff(i)=sum(s2(ismembertol(dsprange,prefdsp,tol))).^beta(1);
        predD60_N_c_sort = (B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^beta(1).*X1(:,1) + ((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^beta(1)).*X1(:,2))./ ...
                (B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^beta(1) + beta(4).*((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^beta(1)) + beta(2)) + beta(3);
            SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_norm(i) = 100*(1 - SSR/SST);
         else
           PV_norm(i) = 0;
         end;
         
       else
       fun = @(w)sum(((B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^w(1).*X1(:,1) +((sum(s2.*normpdf(dsprange,prefdsp,v))).^w(1)).*X1(:,2))./ ...
                      (B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^w(1) + w(4).*((sum(s2.*normpdf(dsprange,prefdsp,v))).^w(1)) + w(2)) + w(3) - Y).^2);
        [beta, func(i)] = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        
        nv(i)=beta(1);
        ncoeff(i)=sum(s1.*normpdf(dsprange,prefdsp,v)).^beta(1);
        fcoeff(i)=sum(s2.*normpdf(dsprange,prefdsp,v)).^beta(1);
        
        predD60_N_c_sort = (B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^beta(1).*X1(:,1) + ((sum(s2.*normpdf(dsprange,prefdsp,v))).^beta(1)).*X1(:,2))./ ...
                (B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^beta(1) + beta(4).*((sum(s2.*normpdf(dsprange,prefdsp,v))).^beta(1)) + beta(2)) + beta(3);
            SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_norm(i) = 100*(1 - SSR/SST);
         else
           PV_norm(i) = 0;
         end;
         
         
         
         
     
%         %add four conditions up
%         for step=lb:s:ub
%        i=i+1;
%        func(i) = sum((((sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1).*X1(:,1) +w(4).*((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)).*X1(:,2))./ ...
%                       ((sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1) + w(4).*((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)) + beta(2)) + beta(3) - Y).^2);
%         
       end
        end
        vvec=(lb1:s:ub1);
        
%         figure(100)
%         clf(100)
%             subplot(1,2,1)
%     plot(vvec,func)
%     hold on
%     title(sprintf(" prefdsp=%f ",prefdsp))
%     subplot(1,2,2)
%     plot(vvec,gradient(func))
%     hold on
%     
%     
%     pause
        
        elseif alphayesno=="YesAF"
            A = [];
        b = [];
        Aeq = []; 
        beq = [];
     beta0 = [1 0 0 0.5]; % beta(1) = n, beta(2) = sigma, beta(3) = C
        
        %look for filter width v within disparity range -1.6 to 1.6 deg
        lb = [0 0 0 0.01];   %alpha 0.1 to10
        ub = [100 5000 100 100];    % changed these upper bounds on 1/9/2020
        lb1 = min(area);   % changed the sigma lower bound from -inf to 0 on 1/9/2020
        ub1 = max(area); %15;    % If we assume 3.2 is the max full width at half max v (or std)=3.2/2.355
        
        s=(max(area)-min(area))/20;
        i=0;
        tol=0.08; %tolerance to prefdsp value incase there is no exact match for prefdsp in dsprange
        %costfun=
        for step=[lb1:s:ub1]
       i=i+1;
         
         tol1=0.01;
         v=min(vec(ismembertol(area,step,tol1)));
       if step==0
           fun = @(w)sum(((B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^w(1).*X1(:,1) +((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^w(1)).*X1(:,2))./ ...
                      (w(4).*B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^w(1) + ((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^w(1)) + w(2)) + w(3) - Y).^2);
        [beta, func(i)] = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        
        nv(i)=beta(1);
        ncoeff(i)=sum(s1(ismembertol(dsprange,prefdsp,tol))).^beta(1);
        fcoeff(i)=sum(s2(ismembertol(dsprange,prefdsp,tol))).^beta(1);
        predD60_N_c_sort = (B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^beta(1).*X1(:,1) + ((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^beta(1)).*X1(:,2))./ ...
                (beta(4).*B*(sum(s1(ismembertol(dsprange,prefdsp,tol)))).^beta(1) + ((sum(s2(ismembertol(dsprange,prefdsp,tol)))).^beta(1)) + beta(2)) + beta(3);
            SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_norm(i) = 100*(1 - SSR/SST);
         else
           PV_norm(i) = 0;
         end;
           
      
       else
       fun = @(w)sum(((B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^w(1).*X1(:,1) +((sum(s2.*normpdf(dsprange,prefdsp,v))).^w(1)).*X1(:,2))./ ...
                      (w(4).*B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^w(1) + ((sum(s2.*normpdf(dsprange,prefdsp,v))).^w(1)) + w(2)) + w(3) - Y).^2);
        [beta, func(i)] = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        
        nv(i)=beta(1);
        ncoeff(i)=sum(s1.*normpdf(dsprange,prefdsp,v)).^beta(1);
        fcoeff(i)=sum(s2.*normpdf(dsprange,prefdsp,v)).^beta(1);
        
        predD60_N_c_sort = (B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^beta(1).*X1(:,1) + ((sum(s2.*normpdf(dsprange,prefdsp,v))).^beta(1)).*X1(:,2))./ ...
                (beta(4).*B*(sum(s1.*normpdf(dsprange,prefdsp,v))).^beta(1) + ((sum(s2.*normpdf(dsprange,prefdsp,v))).^beta(1)) + beta(2)) + beta(3);
            SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_norm(i) = 100*(1 - SSR/SST);
         else
           PV_norm(i) = 0;
         end;
         
         
     
%         %add four conditions up
%         for step=lb:s:ub
%        i=i+1;
%        func(i) = sum((((sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1).*X1(:,1) +w(4).*((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)).*X1(:,2))./ ...
%                       ((sum(s1.*normpdf(dsprange,prefdsp,step))).^beta(1) + w(4).*((sum(s2.*normpdf(dsprange,prefdsp,step))).^beta(1)) + beta(2)) + beta(3) - Y).^2);
%         
       end
        end
         vvec=(lb1:s:ub1);
%         figure(100)
%         clf(100)
%             subplot(1,2,1)
%     plot(vvec,func)
%     hold on
%     title(sprintf(" prefdsp=%f ",prefdsp))
%     subplot(1,2,2)
%     plot(vvec,gradient(func))
%     hold on
%     
%     
%     pause
        
%alpha follows only the nonpreferred dir; split the tuning curves into "Near side" and "far side"      
%not ready fix it
elseif alphayesno=="Yes2"
    beta0 = [1 0 0 0 0.1]; % beta(1) = n, beta(2) = alpha sigma, beta(3) = C beta(4)=alpha beta(5)=w(4)
        A = [];
        b = [];
        Aeq = []; 
        beq = [];
        
        %look for filter width w(4) within disparity range -1 to 1 deg
       lb = [0 0 0 0 0 ];   % changed the sigma lower bound from -5000 to 0 on 1/9/2020         ub = [100 100 5000 100];  
       ub = [100 1 5000 100 3.2];     % changed these upper bounds on 1/9/2020
      fun = @(w)sum((((sum(s1.*normpdf(dsprange,prefdsp,w(5)))).^w(1).*X1(:,1) +((sum(s2.*normpdf(dsprange,prefdsp,w(5)))).^w(1)).*X1(:,2))./ ...          ((sum(s1.*normpdf(dsprange,prefdsp,w(4)))).^w(1) + w(2).*((sum(s2.*normpdf(dsprange,prefdsp,w(4)))).^w(1)) + w(3)) + w(4) - Y).^2);
                      ((sum(s1.*normpdf(dsprange,prefdsp,w(5)))).^w(1) + w(2).*((sum(s2.*normpdf(dsprange,prefdsp,w(5)))).^w(1)) + w(3)) + w(4) - Y).^2);
        beta = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        predD60_N_c_sort = ((sum(s1.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1).*X1(:,1) + ((sum(s2.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1)).*X1(:,2))./ ...
                ((sum(s1.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1) + beta(2).*((sum(s2.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1)) + beta(3)) + beta(4);
        % calculate the weights
         w_normD60_N_c_sort(1) = (sum(s1.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1)./ ...
                ((sum(s1.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1) + beta(2).*((sum(s2.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1)) + beta(3));
         w_normD60_N_c_sort(2) = ((sum(s2.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1))./ ...
                ((sum(s1.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1) + beta(2).*((sum(s2.*normpdf(dsprange,prefdsp,beta(5)))).^beta(1)) + beta(3));
         %HFIT_normD60_N_c_sort(:, whichneuron) = [beta(1)   beta(2)   beta(3)   beta(4) beta(5)]';
        % percentage of w(4)ariance (PV)
         SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_normD60_N_c_sort = 100*(1 - SSR/SST);
         else
           PV_normD60_N_c_sort = 0;
         end;
     
         % R square
        [R,P]=corrcoef(Y, predD60_N_c_sort);
        a = R(2,1)^2;
        Rsq_normD60_N_c_sort = a;
        
elseif alphayesno=="NoYes"
beta0 = [1 0 0 0.1 0.1]; % beta(1) = n, beta(2) = sigma, beta(3) = C
        A = [];
        b = [];
        Aeq = []; 
        beq = [];
        %look for filter width v within disparity range -1.6 to 1.6 deg
        lb = [0 0 0 0 0];   % changed the sigma lower bound from -inf to 0 on 1/9/2020
        ub = [100 5000 100 3.2 10];    % changed these upper bounds on 1/9/2020
       
       fun = @(w)sum(((w(5)*(sum(s1.*normpdf(dsprange,prefdsp,w(4)))).^w(1).*X1(:,1) +((sum(s2.*normpdf(dsprange,prefdsp,w(4)))).^w(1)).*X1(:,2))./ ...
                      (w(5)*(sum(s1.*normpdf(dsprange,prefdsp,w(4)))).^w(1) + ((sum(s2.*normpdf(dsprange,prefdsp,w(4)))).^w(1)) + w(2)) + w(3) - Y).^2);
        beta = fmincon(fun,beta0,A,b,Aeq,beq,lb,ub);
        predD60_N_c_sort = (beta(5)*(sum(s1.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1).*X1(:,1) + ((sum(s2.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1)).*X1(:,2))./ ...
                (beta(5)*(sum(s1.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1) + ((sum(s2.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1)) + beta(2)) + beta(3);
        % calculate the weights
         w_normD60_N_c_sort(1) = beta(5)*(sum(s1.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1)./ ...
                (beta(5)*(sum(s1.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1) + ((sum(s2.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1)) + beta(2));
         w_normD60_N_c_sort(2) = ((sum(s2.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1))./ ...
                (beta(5)*(sum(s1.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1) + ((sum(s2.*normpdf(dsprange,prefdsp,beta(4)))).^beta(1)) + beta(2));
         %HFIT_normD60_N_c_sort(:, whichneuron) = [beta(1)   beta(2)   beta(3) beta(4)]';
        % percentage of w(4)ariance (PV)
         SST = sum((Y-mean(Y)).^2)/length(Y);
         SSR = sum((predD60_N_c_sort-Y).^2)/length(predD60_N_c_sort);
         if SSR<=SST
           PV_normD60_N_c_sort = 100*(1 - SSR/SST);
         else
           PV_normD60_N_c_sort = 0;
         end;
     
         % R square
        [R,P]=corrcoef(Y, predD60_N_c_sort);
        a = R(2,1)^2;
        Rsq_normD60_N_c_sort = a;

end



