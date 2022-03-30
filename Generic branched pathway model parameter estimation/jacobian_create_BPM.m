clc
clear all



%%% METABOLITE SYMBOLS

syms x0 x1 x2 x3 x4 g1 g2 g3 g4 g5 g6 f13 f21 f32 f43 f44 f51 f64;


xdot(1) = (g1*x0)/(x3^f13) - g2*x1^f21 - g5*x1^f51;
xdot(2) = g2*x1^f21 - g3*x2^f32;
xdot(3) = g3*x2^f32 - g4*x3^f43*x4^f44;
xdot(4) = g5*x1^f51 - g6*x4^f64;
xdot(5) = 0;
xdot(6) = 0;
xdot(7) = 0;
xdot(8) = 0;
xdot(9) = 0;
xdot(10) = 0;
xdot(11) = 0;




var = [x1, x2, x3, x4, f13, f21, f32, f43, f44, f51, f64];


C = {'x1', 'x2', 'x3', 'x4', 'f13', 'f21', 'f32', 'f43', 'f44', 'f51', 'f64'};

diary BPM_parameter_estimation_jac.m
diary off


str =  ["(g1*x0)/(x3^f13) - g2*x1^f21 - g5*x1^f51";
"g2*x1^f21 - g3*x2^f32";
"g3*x2^f32 - g4*x3^f43*x4^f44";
"g5*x1^f51 - g6*x4^f64";
"0";
"0";
"0";
"0";
"0";
"0";
"0"];



count = 0;
for i=1:11
   expr = xdot(i); 
   for j=1:11
      s = contains(str(i),C(j));  
      if(s==1)         
         exp_diff = diff(expr, var(j));
         if(exp_diff ~= 0)
            simp_expr = simplify(exp_diff); 
            diary on
            fprintf("\njac(%d, %d) = %s",i,j,simp_expr);
%             disp(simp_expr); 
            disp(";");
            diary off
            count = count + 1;
         end
%       else
%          diary on
%          fprintf("\n\njac(%d, %d) = ",i,j);
%          disp('0');
%          diary off 
      end
   end
end
   
diary on
fprintf("\n\n################\n\n Number of non zero element = %d",count);
diary off
   
   









% % % %%% METABOLITE SYMBOLS
% % % 
% % % syms x0 x1 x2 x3 x4 g1 g2 g3 g4 g5 g6 f13 f21 f32 f43 f44 f51 f64 g1 g2 g3 g4 g5 g6;
% % % 
% % % 
% % % xdot(1) = (g1*x0)/(x3^f13) - g2*x1^f21 - g5*x1^f51;
% % % xdot(2) = g2*x1^f21 - g3*x2^f32;
% % % xdot(3) = g3*x2^f32 - g4*x3^f43*x4^f44;
% % % xdot(4) = g5*x1^f51 - g6*x4^f64;
% % % xdot(5) = 0;
% % % xdot(6) = 0;
% % % xdot(7) = 0;
% % % xdot(8) = 0;
% % % xdot(9) = 0;
% % % xdot(10) = 0;
% % % xdot(11) = 0;
% % % xdot(12) = 0;
% % % xdot(13) = 0;
% % % xdot(14) = 0;
% % % xdot(15) = 0;
% % % xdot(16) = 0;
% % % xdot(17) = 0;
% % % 
% % % 
% % % 
% % % 
% % % var = [x1, x2, x3, x4, f13, f21, f32, f43, f44, f51, f64, g1, g2, g3, g4, g5, g6];
% % % 
% % % 
% % % C = {'x1', 'x2', 'x3', 'x4', 'f13', 'f21', 'f32', 'f43', 'f44', 'f51', 'f64', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6'};
% % % 
% % % diary BPM_parameter_estimation_jac.m
% % % diary off
% % % 
% % % 
% % % str =  ["(g1*x0)/(x3^f13) - g2*x1^f21 - g5*x1^f51";
% % % "g2*x1^f21 - g3*x2^f32";
% % % "g3*x2^f32 - g4*x3^f43*x4^f44";
% % % "g5*x1^f51 - g6*x4^f64";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0";
% % % "0"];
% % % 
% % % 
% % % 
% % % count = 0;
% % % for i=1:17
% % %    expr = xdot(i); 
% % %    for j=1:17
% % %       s = contains(str(i),C(j));  
% % %       if(s==1)         
% % %          exp_diff = diff(expr, var(j));
% % %          if(exp_diff ~= 0)
% % %             simp_expr = simplify(exp_diff); 
% % %             diary on
% % %             fprintf("\njac(%d, %d) = %s",i,j,simp_expr);
% % % %             disp(simp_expr); 
% % %             disp(";");
% % %             diary off
% % %             count = count + 1;
% % %          end
% % % %       else
% % % %          diary on
% % % %          fprintf("\n\njac(%d, %d) = ",i,j);
% % % %          disp('0');
% % % %          diary off 
% % %       end
% % %    end
% % % end
% % %    
% % % diary on
% % % fprintf("\n\n################\n\n Number of non zero element = %d",count);
% % % diary off






   
