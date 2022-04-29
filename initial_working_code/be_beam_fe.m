function [ mass_assem, stif_assem] = be_beam_fe( n, m, ei, l, alpha )
me = m/n;
le = l/n;
[ mass_matr, stif_matr] = be_beam_matr( me, ei, le, alpha );
mass_assem_temp = zeros(2 * n + 2, 2 * n + 2);
stif_assem_temp = zeros(2 * n + 2, 2 * n + 2);
for j = 1:n
mass_assem_temp(2 * j - 1 : 2 * j + 2, 2 * j - 1 : 2 * j + 2) = mass_assem_temp(2 * j - 1 : 2 * j + 2, 2 * j - 1 : 2 * j + 2) + mass_matr; % add overlapping block diagonal matrices
stif_assem_temp(2 * j - 1 : 2 * j + 2, 2 * j - 1 : 2 * j + 2) = stif_assem_temp(2 * j - 1 : 2 * j + 2, 2 * j - 1 : 2 * j + 2) + stif_matr; % add overlapping block diagonal matrices
end

mass_assem = mass_assem_temp([2:end-2,size(mass_assem_temp,1)],[2:end-2,size(mass_assem_temp,2)]); % remove 1st row/col, second to last row/col
stif_assem = stif_assem_temp([2:end-2,size(stif_assem_temp,1)],[2:end-2,size(stif_assem_temp,2)]); % remove 1st row/col, second to last row/col
%   1 2 3 4 5 6 7 8
%1 |       |
%2 |       |
%3 |   |   |   |
%4 |   |   |   |
%5     |       |
%6     |       |
%
%
%
