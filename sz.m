function x = sz(lambda_0, delta_d_0, lambda_1, c_t, f_0_t, T, x_t) 
   if ~exist('x_t','var')
      x_t = 0:1:T;
   end

    lambda_0_c = @(c) lambda_0 - delta_d_0.*c./(c+1);

    function y = fin_size(t)
        c = interp1(x_t,c_t,t);
        f_0 = interp1(x_t,f_0_t,t);

        y = (lambda_0_c(c)-lambda_1).*f_0+(lambda_1);
    end

    x = integral(@fin_size,min(x_t),T);
end
