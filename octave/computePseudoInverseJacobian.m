function Jinv = computePseudoInverseJacobian(Jnum, rankThreshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Jinv = computePseudoInverseJacobian(Jnum)
% Task: Compute the pseudo inverse of the Jacobian matrix and correct it if loosing rank
%
% Inputs:
%	- Jnum: x-coordinate (in m)
%	- rankThreshold: threshold value to estimate if the matrix is losing rank
%
% Outputs: 
%	- Jinv: pseudo inverse of the Jacobian matrix
%	
%
% author: Guillaume Gibert, guillaume.gibert@ecam.fr
% date: 14/09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[U, S, V] = svd(Jnum);

% check for rank and replace singular values lower than threshold by 0
S_inv = S';
for l_singVal=size(Jnum,2):-1:1
	if (S(l_singVal,l_singVal)/S(1,1) < rankThreshold)
		S_inv(l_singVal,l_singVal) = 0;
		%disp('[INFO](computePseudoInverseJacobian) Jacobian matrix loosing rank!')
	else
		S_inv(l_singVal,l_singVal) = 1.0/S(l_singVal,l_singVal);
	end
end

Jinv = V*S_inv*U';


