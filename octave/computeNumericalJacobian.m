function Jnum = computeNumericalJacobian(theta, d, a, alpha, revolute)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Jnum = computeNumericalJacobian(theta, d, a, alpha)
% Task: Compute Jnum (i.e. numerical Jacobian matrix of the robot) given the D-H parameters
%
% Inputs:
%	- theta: an array of theta parameters (rotation around z in degrees)
%	- d: an array of d parameters (translation along z in mm)
%	- a: an array of a parameters (translation along x in mm)
%	- alpha: an array of alpha parameters (rotation around x in degrees)
%	- revolute: an array of value equal to 1 if the joint is revolute and 0 if prismatic
%
% Outputs: 
%	- Jnum: numerical Jacobian matrix of the robot
%	
%
% author: Guillaume Gibert, guillaume.gibert@ecam.fr
% date: 14/09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jnum = [];

for l_joint=1:size(theta,1)
	jTee = dh2ForwardKinematics(theta, d, a, alpha, l_joint);
	jtee = jTee(1:3, end);
	eeTj = inverse3DTransformationMatrix(jTee);
	eeRj = eeTj(1:3,1:3);
	% compute Tj 
	Tj = [eeRj -eeRj*skewSymetricMatrix(jtee(1), jtee(2), jtee(3));
		zeros(3,3) eeRj];
		
	if (revolute(l_joint)==1)
		Jnum = [Jnum Tj(:,6)]; % revolute joint
	elseif (revolute(l_joint)==2)
		Jnum = [Jnum Tj(:,3)]; % prismatic joint
	end
end
