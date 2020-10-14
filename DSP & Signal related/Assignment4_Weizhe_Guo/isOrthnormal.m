function [isOrth] = isOrthnormal(B)
   prod = B'*B;
   if max(max(abs(prod - eye(size(prod)))))<10000*eps
       isOrth = 1;
   else
       isOrth = 0;
   end
end