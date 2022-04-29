function [ alfa, beta] = rayparam(omega, zeta)
 % param = lsqr(transpose([1 / (2 * omega), omega / 2]), zeta * ones(size(omega))); % ne treba transpose?
 [param,~] = lsqr(([1 ./ (2 * omega), omega / 2]), zeta * ones(size(omega)));
 alfa = real(param(1));
 beta = real(param(2));

end