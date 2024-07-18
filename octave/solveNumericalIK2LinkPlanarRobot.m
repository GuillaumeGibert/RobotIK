function solveNumericalIK2LinkPlanarRobot(L1, L2, q1, q2, x_target, y_target, pController, rankThreshold, distanceThreshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function solveNumericalIK2LinkPlanarRobot(L1, L2, q1, q2, x_target, y_target, pController, rankThreshold, distanceThreshold)
% ex.: solveNumericalIK2LinkPlanarRobot(2.0, 1.0, 100, 50, -1.5, 0.1, 0.5, 0.2, 0.01)
% Task: solve Inverse Kinematics (if it exists) numerically for a 2-link planar robot and compare it to the analitical one
%
% Inputs:
%	- L1: length of link 1 (in m)
%	- L2: length of link 1 (in m)
%	- q1: value of the first joint angle (in degrees)
%	- q2:  value of the second joint angle (in degrees)
%	- x_target: target x-coordinate (in m)
%	- y_target:  target y-coordinate (in m)
%	- pController: the value of the proportionnal controller when computing Xdot (a.u.)
%	- rankThreshold: threshold value to determine if the pseudoinverse matrix is close to losing rank (a.u.)
%	- distanceThreshold:   threshold value to stop the optimization loop (in m)
%
% Output: None
%	
%
% author: Guillaume Gibert, guillaume.gibert@ecam.fr
% date: 09/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creates the D-H parameters from the 2-link planar robot settings
theta = [q1 q2]';
d = [0 0]';
a = [L1 L2]';
alpha = [0 0]';
revolute = [1 1];

% creates the target point in Cartesian space in homogeneous format
X_target = [x_target y_target 0]';

% calls the analytical method
fprintf("[Analytical method] Launching the analytical method...");
[nbSol, qi] = solveAnalyticalIK2LinkPlanarRobot(L1, L2, X_target(1), X_target(2));
fprintf(" Done!\n");

if (nbSol == 0)
	fprintf("[Analytical method] There is no solution!\n");
	return;
elseif (nbSol == 1)
	fprintf("[Analytical method] There is one solution: (q1, q2) = (%f, %f)\n", qi(1,1), qi(1,2));
elseif (nbSol == 2)
	fprintf("[Analytical method] There are two solutions: \n");
	fprintf("(q1, q2) = (%f, %f)\n", qi(1,1), qi(1,2));
	fprintf("(q1, q2) = (%f, %f)\n", qi(2,1), qi(2,2));
end


% calls the numerical method

% creates a "boolean" variable to control the optimization loop
estimationInProgress = 1;
nbIterations = 0;

% creates a figure
figure;

% starts the optimization loop
fprintf("[Numerical method] Launching the numerical method...");
while (estimationInProgress)

	%1- Computes Vee
	Vee = computeVee(X_target, pController, theta, d, a, alpha);
	%2- Computes Numerical Jacobian
	Jnum = computeNumericalJacobian(theta, d, a, alpha,revolute); 
	%3- compute Pseudo Inverse of Jacobian
	Jinv = computePseudoInverseJacobian(Jnum, rankThreshold);
	%4- Estimates deltaQ
	deltaQ = Jinv * Vee;
	%5- Updates qi values
	for l_joint=1:size(theta,1)
		theta(l_joint) +=  deltaQ(l_joint)/pController;
	end
	%6- Updates the current position of the ee
	jTee = dh2ForwardKinematics(theta, d, a, alpha, 1);
	X_current = jTee*[0 0 0 1]';
	X_current(end) = [];
	%7- Estimates the distance with the target position
	distance = 0;
	for l_coord=1:3
		distance += (X_target(l_coord)-X_current(l_coord))^2;
	end
	distance = sqrt(distance);
	%8- Checks if the loop must still run 
	if (distance < distanceThreshold)
		estimationInProgress = 0;
	end
	%9- Estimates the position of the tip of J2 = ee
	wTee = dh2ForwardKinematics(theta, d, a, alpha,1);
	wPee = wTee * [0 0 0 1]';
		
	%fprintf("[Numerical method] Iteration #%d, Distance: %f m\n", nbIterations, distance);
	
	if (1)
	subplot(1,3,1); % left panel
		drawCircle(0, 0, L1+L2); % draw the outer circle of the workspace
		drawCircle(0, 0, L1-L2); % draw the inner circle of the workspace
		xlabel('x-coord (m)');
		ylabel('y-coord (m)');
		plot(x_target,y_target,'or'); hold on; % target
		plot(wPee(1,1),wPee(2,1),'*'); hold on;
		xlim([-(L1+L2) L1+L2]);
		ylim([-(L1+L2) L1+L2]);
	subplot(1,3,2); % middle panel
		xlabel('Iteration (a.u.)');
		ylabel('qi value (°)');
		plot(nbIterations, theta(1), 'b'); hold on;
		plot(nbIterations, theta(2), 'r'); hold on;
		%xlim([0 2000]);
		ylim([-180 180]);
	subplot(1,3,3); % middle panel
		xlabel('Iteration (a.u.)');
		ylabel('Distance to target (m)');
		plot(nbIterations, distance, 'g'); hold on;
		plot([1 nbIterations], [distanceThreshold distanceThreshold], 'r'); hold on;
		%xlim([0 2000]);
		ylim([0 2*(L1+L2)]);
	end
	
	nbIterations++;
	pause(0.001);
	
end

fprintf(" Done!\n");
fprintf("[Numerical method] The solution is: (q1, q2) = (%f, %f)\n", theta(1), theta(2));
