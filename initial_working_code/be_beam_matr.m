function [ mass_matr, stif_matr] = be_beam_matr( m, ei, l, alpha )
mass_matr = m*diag([0.5, alpha*l^2, 0.5, alpha*l^2]);
stif_matr = ei / l^3 * [ 12,   6*l,  -12,   6*l;
	                    6*l, 4*l^2, -6*l, 2*l^2;
						-12,  -6*l,   12,  -6*l;
						6*l, 2*l^2, -6*l, 4*l^2];

end