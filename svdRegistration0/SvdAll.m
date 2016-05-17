function [Uall,Un,Sv] = SvdAll(Uall)

Uall = reshape(Uall,[],size(Uall,3)*size(Uall,4));
[Un, Sv, ~] = svdecon(Uall'*Uall);
Sv = diag(Sv);
