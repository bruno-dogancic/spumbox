function [A, B, C, D] = mehss(M, P, K, B1, C1)
[Adss, Bdss, C, D, E] = mehdss(M, P, K, B1, C1);
A = linsolve(E, Adss);
B = linsolve(E, Bdss);