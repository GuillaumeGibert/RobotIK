function ssm = skewSymetricMatrix(x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ssm = skewSymetricMatrix(x,y,z)
% Task: Create a skew symetric matrix from a 3D vector coordinates
%
% Inputs:
%	- x: x-coordinate (in m)
%	- y: y-coordinate  (in m)
%	- z: z-coordinate  (in m)
%
% Outputs: 
%	- ssm: the skew symetric matrix of the input vector
%	
%
% author: Guillaume Gibert, guillaume.gibert@ecam.fr
% date: 14/09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ssm = [	0 -z y;
		z 0 -x;
		-y x 0];