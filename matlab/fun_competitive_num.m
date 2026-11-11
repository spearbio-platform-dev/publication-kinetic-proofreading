function vfE = fun_competitive_num(vfCA,PARAs)
%FUN_COMPETITIVE_NUMSIM Summary of this function goes here
%   Detailed explanation goes here
vfE = nan(size(vfCA));
CB = PARAs(1);
CE = PARAs(2);
K_D = PARAs(3);
K_DE = PARAs(4);

syms A B E BA BE

for i = 1:length(vfCA)
    CA = vfCA(i);
    eqns = [ ...
        BA/(A*B) == 1/K_D, ...
        BE/(B*E) == 1/K_DE, ...
        A + BA   == CA, ...
        B + BE + BA == CB, ...
        E + BE   == CE];
    S = vpasolve(eqns, [A B E BA BE], [CA CB CE 0 0]);
    vfE(i) = max(double(S.E));
end





end

