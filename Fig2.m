   lambda_0 = 0.04; delta_d_0 = 0.08;
   k_0 = 0.0004; m_0 = 0.0004;
   mu = 0.0004; k = [k_0,0,k_0];
   lambda_1 = 0.001;
   nu = 0.004; m = [0,m_0,m_0];

   max_dose = 10; 
   T = 1200;
   x_t = 0:1:T;
   difference = zeros(3,size(x_t,2));

   figure;
   tiledlayout(4,3);

   for ell=1:3
       nexttile(ell);

       ic_in = get_equilibrium(0,lambda_0,delta_d_0,lambda_1,mu,k(ell),nu,m(ell));
       [a,b,c,d,i,j,e,f,g] = optimal_control_baedi_linur(lambda_0,delta_d_0,lambda_1,mu,k(ell),nu,m(ell),T,max_dose,ic_in);
    
       results = oc_maeling;
    
       results.c_t = a;
       results.c_t_iter = b;
       results.c_t_orig_iter = c;
       results.f_0_t = d;
       results.f_0_t_iter = i;
       results.f_0_t_orig_iter = j;
       results.gamma_t = e;
       results.final_size = f;
       results.final_size_vec = g;
    
       plot([0:1:T],results.c_t,'Color',[78 107 166]/255,'LineWidth',3)
       set(gca,'Yscale','log')
       ylim([10^(-0) max_dose])
    
       [~, c_l] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k(ell), nu, m(ell), max_dose);
       hold on
       yline(c_l,'--black','LineWidth',3);

       color = 0.3;
       rectangle('Position',[900,1,300,10],'FaceColor',color*ones(1,3),'EdgeColor',color*ones(1,3),'LineWidth',1,'FaceAlpha',.9)

       set(gca,'fontsize', 14);
       xlabel('Time $t$','Interpreter','Latex','FontSize',19);
       ylabel('Dose $c$ (relative to $EC_{50}^d$)','Interpreter','Latex','FontSize',19);
       legend('Optimal treatment','Optimal equilibrium','Interpreter','Latex','FontSize',19);
       set(gca,'xlim',[0 T]);

       nexttile(ell+3);
       f_0 = get_equilibrium(c_l,lambda_0,delta_d_0,lambda_1,mu,k(ell),nu,m(ell));
       plot(0:1:T,results.f_0_t,'Color',[212 10 0]/255,'LineWidth',3);
       hold on
       yline(f_0,'--black','LineWidth',3);

       color = 0.3;
       rectangle('Position',[900,0,300,1],'FaceColor',color*ones(1,3),'EdgeColor',color*ones(1,3),'LineWidth',1,'FaceAlpha',.9)

       set(gca,'fontsize', 14);
       xlabel('Time $t$','Interpreter','Latex','FontSize',19);
       ylabel('Sensitive cell proportion $f_0$','Interpreter','Latex','FontSize',19);
       legend('Optimal treatment','Optimal equilibrium','Interpreter','Latex','FontSize',19);
       set(gca,'xlim',[0 T]);
       ylim([0 1]);

       nexttile(ell+6);
       vec = zeros(1,size(x_t,2));
       for err=1:size(vec,2)
          vec(err) = sz(lambda_0, delta_d_0, lambda_1, results.c_t, results.f_0_t, x_t(err), x_t);
       end
       vec = 10^10*exp(vec);
       plot(x_t,vec,'Color',[78 107 166]/255,'LineWidth',3);
       set(gca,'Yscale','log');
       hold on
       
       lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
       mu_c = @(c) mu + k(ell).*c;
       nu_c = @(c) nu - m(ell).*c;
       [rho_l, c_l] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k(ell), nu, m(ell), max_dose);
       c = c_l;
       A = [lambda_0_c(c)-mu_c(c), mu_c(c); nu_c(c),lambda_1-nu_c(c)];
       vec_2 = zeros(1,size(x_t,2));
       for err=1:size(vec_2,2)
           vec_2(err) = 10^10*sum([ic_in,1-ic_in]*expm(x_t(err)*A));
           difference(ell,err) = 1-vec(err)/vec_2(err);
       end
       plot(x_t,vec_2,'Color',[212 10 0]/255,'LineWidth',3);
       xlim([0 T]);
       color = 1;

       set(gca,'fontsize', 14);
       xlabel('Time $t$','Interpreter','Latex','FontSize',19);
       ylabel('Tumor size $n(t)$ (log scale)','Interpreter','Latex','FontSize',19);
       legend('Optimal treatment','Constant dose','Interpreter','Latex','FontSize',19)
       ylim([0.4*10^8 10^10])
       rectangle('Position',[600,0,400,10^10-10],'FaceColor',color*ones(1,3),'EdgeColor',color*ones(1,3),'LineWidth',1,'FaceAlpha',1);
   end

    ic_in = [1,0.8,0.6,0.4,0.2,0.05];

   for ell=1:3
       nexttile;
       hold on
       for err = 1:size(ic_in,2)
           [a,b,c,d,i,j,e,f,g] = optimal_control_baedi_linur(lambda_0,delta_d_0,lambda_1,mu,k(ell),nu,m(ell),T,max_dose,ic_in(err));
        
           results = oc_maeling;
        
           results.c_t = a;
           results.c_t_iter = b;
           results.c_t_orig_iter = c;
           results.f_0_t = d;
           results.f_0_t_iter = i;
           results.f_0_t_orig_iter = j;
           results.gamma_t = e;
           results.final_size = f;
           results.final_size_vec = g;
        
           plot([0:1:T],results.c_t,'LineWidth',3)
           set(gca,'Yscale','log')
           ylim([10^(-1) max_dose])
    
           color = 0.3;
           rectangle('Position',[900,1,300,99],'FaceColor',color*ones(1,3),'EdgeColor',color*ones(1,3),'LineWidth',1,'FaceAlpha',.9)
    
           set(gca,'fontsize', 14);
           xlabel('Time $t$','Interpreter','Latex','FontSize',19);
           ylabel('Dose $c$ (relative to $EC_{50}^d$)','Interpreter','Latex','FontSize',19);
           set(gca,'xlim',[0 T]);
       end
       [~, c_l] = best_constant_dose_limits(max_dose, lambda_0, delta_d_0, lambda_1, mu, k(ell), nu, m(ell), max_dose);
       hold on
       yline(c_l,'--black','LineWidth',3);
       legend({'$f_0(0)=1$','$f_0(0)=0.8$','$f_0(0)=0.6$','$f_0(0)=0.4$','$f_0(0)=0.2$','$f_0(0)=0.01$',''},'Location','northeast','NumColumns',2,'Interpreter','Latex','FontSize',19);
       xlim([0 600]);
   end

   function f_0 = get_equilibrium(c,lambda_0,delta_d_0,lambda_1,mu,k,nu,m) 
      lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);
      mu_c = @(c) mu + k.*c;
      nu_c = @(c) nu - m.*c;

      A = [lambda_0_c(c)-mu_c(c), mu_c(c); nu_c(c),lambda_1-nu_c(c)];
      
      sigma = max(eig(A));    
      f_0 = A(2,1)/(sigma-A(1,1)+A(2,1));
   end