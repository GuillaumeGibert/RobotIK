function Vee = computeVee(X_target, pController, theta, d, a, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Vee = computeVee(X_target, theta, d, a, alpha)
% Task: Compute Vee (i.e. velocity of the end-effector in {ee})
%
% Inputs:
%	- X_target: target position/orientation vector of the end-effector in {ee} (in mm and degrees)
%	- pController: p value of the proportional controller (a.u.)
%	- theta: an array of theta parameters (rotation around z in degrees)
%	- d: an array of d parameters (translation along z in mm)
%	- a: an array of a parameters (translation along x in mm)
%	- alpha: an array of alpha parameters (rotation around x in degrees)
%
% Outputs: 
%	- Vee: Velocity of the end-effector in {ee}
%	
%
% author: Guillaume Gibert, guillaume.gibert@ecam.fr
% date: 14/09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%0th step: current X
jTee = dh2ForwardKinematics(theta, d, a, alpha, 1);
X_current = jTee*[0 0 0 1]';
X_current(end) = [];

% 1st step: deltaX
if (size(X_target,1) ~= size(X_current,1))
	disp('[ERROR](computeVee)-> sizes of input arrays do not match!')
	return;
end
deltaX = X_target-X_current;
deltaX = [deltaX; 0; 0; 0]; % add zeros for orientation as we want only to control position

% 2nd step: Xdot
Xdot = pController*deltaX;

% 3rd step: eeRb
jTee = dh2ForwardKinematics(theta, d, a, alpha, 1);
eeRb = jTee(1:3, 1:3)';

% 4th step: Vee
Vee = [	eeRb zeros(3,3);
		zeros(3,3) eeRb] * Xdot;
		