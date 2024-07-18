#include "Kinematics.h"

Kinematics::Kinematics():
	m_bAreDHParametersSet(false), m_bIsDatabaseLoaded(false),
	m_i32NbJoints(0), 
	m_i32MaxNbIterations(0), m_f64PController(0.0), m_f64RankThreshold(0.0), m_f64DistanceThreshold(0.0)
{}

Kinematics::Kinematics(std::vector<double>& vDHTheta, std::vector<double>& vDHD, std::vector<double>& vDHA, std::vector<double>& vDHAlpha, std::vector<int>& vJointType, std::vector<double>& vJointOffset, std::vector<double>& vQiValues, int i32NbJoints, int i32MaxNbIterations, double f64PController, double f64RankThreshold, double f64DistanceThreshold):
	m_bAreDHParametersSet(true), m_bIsDatabaseLoaded(false),
	m_vDHTheta(vDHTheta), m_vDHD(vDHD), m_vDHA(vDHA), m_vDHAlpha(vDHAlpha), m_vJointType(vJointType), m_vJointOffset(vJointOffset), m_vQiValues(vQiValues), m_vQiOriginalValues(vQiValues),
	m_i32NbJoints(i32NbJoints), 
	m_i32MaxNbIterations(i32MaxNbIterations), m_f64PController(f64PController), m_f64RankThreshold(f64RankThreshold), m_f64DistanceThreshold(f64DistanceThreshold)
{
	updateDHParameters();
}

Kinematics::~Kinematics()
{
	m_oRandomJointValues.release();
	m_oRandomEndEffectorPositions.release();
}

void Kinematics::updateDHParameters(std::vector<double> vQiValues)
{
	m_vQiValues = vQiValues;
	updateDHParameters();
}

void Kinematics::updateDHParameters()
{
	int l_i32JointIndex = 0;
	
	for (int l_joint = 0; l_joint < m_vJointType.size(); l_joint++)
	{
		if (m_vJointType[l_joint] == KINEMATICS_REVOLUTE_JOINT)
		{
			/*std::cout << "l_i32JointIndex = " << l_i32JointIndex << std::endl;
			std::cout << "l_joint = " << l_joint << std::endl;
			std::cout << "m_vQiValues[l_i32JointIndex] = " << m_vQiValues[l_i32JointIndex] << std::endl;
			std::cout << "m_vJointOffset[l_joint]= " << m_vJointOffset[l_joint] << std::endl;
			*/
			m_vDHTheta[l_joint] = m_vQiValues[l_i32JointIndex] + m_vJointOffset[l_joint];
			l_i32JointIndex++;
		}
		else if (m_vJointType[l_joint] == KINEMATICS_PRISMATIC_JOINT)
		{
			m_vDHD[l_joint] = m_vQiValues[l_i32JointIndex] + m_vJointOffset[l_joint];
			l_i32JointIndex++;
		}
		else //TODO find a way to avoid this else
		{
			m_vDHTheta[l_joint] = m_vJointOffset[l_joint];
		}
	}
	
	/*for (int l_joint = 0; l_joint < m_vDHTheta.size(); l_joint ++)
	{
		std::cout << "(" << m_vDHTheta[l_joint] << ", " << m_vDHD[l_joint] << ", " << m_vDHA[l_joint] << ", " << m_vDHAlpha[l_joint] << ") " << std::endl; 
	}*/
	
}

void Kinematics::setDHParameters(std::vector<double>& vDHTheta, std::vector<double>& vDHD, std::vector<double>& vDHA, std::vector<double>& vDHAlpha, std::vector<int>& vJointType, std::vector<double>& vJointOffset, std::vector<double>& vQiValues, int i32NbJoints)
{
	m_vDHTheta = vDHTheta;
	m_vDHD = vDHD;
	m_vDHA = vDHA;
	m_vDHAlpha = vDHAlpha;
	m_vJointType = vJointType;
	m_vJointOffset = vJointOffset;
	m_vQiValues = vQiValues;
	m_i32NbJoints = i32NbJoints;

	// check if the vectors have the same size
	if (m_vDHTheta.size() != m_vDHD.size() || m_vDHTheta.size() != m_vDHA.size() || m_vDHTheta.size() != m_vDHAlpha.size() || m_vDHTheta.size() != m_vJointType.size() || m_vDHTheta.size() != m_vJointOffset.size())
	{
		std::cerr << "(Kinematics::setDHParameters) The vectors of DH parameters do not have the same size!" << std::endl;
		m_bAreDHParametersSet = false;
	}
	else
	{
		updateDHParameters();
		m_bAreDHParametersSet = true;
	}
}

std::vector<double> Kinematics::getQiValues()
{
	return m_vQiValues;
}

cv::Mat Kinematics::dh2FK(int i32StartJointNumber)
{
	cv::Mat l_oTransformationMatrix_jTee = cv::Mat::eye(4, 4, CV_64F);

	if (!m_bAreDHParametersSet)
	{
		std::cerr << "(Kinematics::dh2FK) D-H parameters have not been set yet!" << std::endl;
		return l_oTransformationMatrix_jTee;
	}

	if (i32StartJointNumber < 0 || i32StartJointNumber > m_vDHTheta.size())
	{
		std::cerr << "(Kinematics::dh2FK) the value of i32StartJointNumber is not in the right range!, i32StartJointNumber= "
			<< i32StartJointNumber << ", [0, " << m_vDHTheta.size() << "[" << std::endl;
		return l_oTransformationMatrix_jTee;
	}


	for (int l_joint = i32StartJointNumber; l_joint < m_vDHTheta.size(); l_joint++)
	{
		TransformationMatrix l_oThetaTransformationMatrix(0.0f, 0.0f, m_vDHTheta[l_joint], 1, 0.0f, 0.0f, 0.0f); //Rz
		TransformationMatrix l_oDTransformationMatrix(0.0f, 0.0f, 0.0f, 1, 0.0f, 0.0f, m_vDHD[l_joint]); //Tz
		TransformationMatrix l_oATransformationMatrix(0.0f, 0.0f, 0.0f, 1, m_vDHA[l_joint], 0.0f, 0.0f); //Tx
		TransformationMatrix l_oAlphaTransformationMatrix(m_vDHAlpha[l_joint], 0.0f, 0.0f, 1, 0.0f, 0.0f, 0.0f); //Rx

		l_oTransformationMatrix_jTee *= l_oThetaTransformationMatrix.getTransformationMat() * l_oDTransformationMatrix.getTransformationMat() * l_oATransformationMatrix.getTransformationMat() * l_oAlphaTransformationMatrix.getTransformationMat();
	}

	return l_oTransformationMatrix_jTee;
}


cv::Mat Kinematics::skewSymmetricMatrix(double fX, double fY, double fZ)
{
	cv::Mat l_oSkewSymmetricMatrix(3, 3, CV_64F, 0.0f);

	l_oSkewSymmetricMatrix.at<double>(0, 0) = 0.0f; l_oSkewSymmetricMatrix.at<double>(0, 1) = -fZ; l_oSkewSymmetricMatrix.at<double>(0, 2) = fY;
	l_oSkewSymmetricMatrix.at<double>(1, 0) = fZ; l_oSkewSymmetricMatrix.at<double>(1, 1) = 0.0f; l_oSkewSymmetricMatrix.at<double>(1, 2) = -fX;
	l_oSkewSymmetricMatrix.at<double>(2, 0) = -fY; l_oSkewSymmetricMatrix.at<double>(2, 1) = fX; l_oSkewSymmetricMatrix.at<double>(2, 2) = 0.0f;

	return l_oSkewSymmetricMatrix;
}

cv::Mat Kinematics::computeVee(cv::Mat oXtarget)
{
	cv::Mat l_oVee(6, 1, CV_64F, 0.0f);
	
	// estimates the current X (position and orientation of the end-effector)
	cv::Mat l_oXcurrent_hc = computeCurrentEndEffectorPosition();

	cv::Mat l_oXcurrent(3, 1, CV_64F, 0.0f);
	l_oXcurrent_hc(cv::Rect(0, 0, l_oXcurrent_hc.cols, l_oXcurrent_hc.rows-1)).copyTo(l_oXcurrent); // remove homogeneous coordinate part i.e. last element (= 1)
	
	// checks if the matrices have the same size
	if (l_oXcurrent.rows != oXtarget.rows || l_oXcurrent.cols != oXtarget.cols)
	{
		std::cerr << "(Kinematics::computeVee) Xcurrent and Xtarget matrices do not have the same size!" << std::endl;
		return l_oVee;
	}

	// computes deltaX
	cv::Mat l_oDeltaX = oXtarget - l_oXcurrent;

	// computes Xdot
	cv::Mat l_oXdot(6, 1, CV_64F, 0.0f); // adds zeros for orientation as we are controlling position only
	l_oDeltaX.copyTo(l_oXdot(cv::Rect(0, 0, l_oDeltaX.cols, l_oDeltaX.rows)));
	l_oXdot *= m_f64PController;

	// retrieves eeRb from bTee
	cv::Mat l_oTransformationMatrix_bTee = dh2FK(0);
	TransformationMatrix l_oTransformationMatrix(l_oTransformationMatrix_bTee);
	cv::Mat l_oeeTb = l_oTransformationMatrix.getTransformationMatInverse(); 
	cv::Mat l_oeeRb = l_oeeTb(cv::Rect(0, 0, 3, 3));

	// creates the transformation matrix Xdot->Vee i.e. [eeRb(3x3) zeros(3x3); zeros(3x3) eeRb(3x3)]
	cv::Mat l_oXdot2Vee(6, 6, CV_64F, 0.0f);
	l_oeeRb.copyTo(l_oXdot2Vee(cv::Rect(0, 0, l_oeeRb.cols, l_oeeRb.rows)));
	l_oeeRb.copyTo(l_oXdot2Vee(cv::Rect(3, 3, l_oeeRb.cols, l_oeeRb.rows)));

	// computes Vee = Xdot2Vee * Xdot
	l_oVee = l_oXdot2Vee * l_oXdot;
		
	return l_oVee;
}

cv::Mat Kinematics::computeNumericalJacobian()
{
	cv::Mat l_oNumericalJacobian(6, m_i32NbJoints, CV_64F, 0.0f);

	int l_i32JointIndex = 0;
	for (int l_joint = 0; l_joint < m_vJointType.size(); l_joint++)
	{
		// checks if the current index corresponds to either a revolute or prismatic joint
		if ((m_vJointType[l_joint] == KINEMATICS_REVOLUTE_JOINT) || (m_vJointType[l_joint] == KINEMATICS_PRISMATIC_JOINT))
		{
			// computes the transformation matrix jTee from the DH parameters
			cv::Mat l_oTransformationMatrix_jTee = dh2FK(l_joint);
			TransformationMatrix l_oTransformationMatrix(l_oTransformationMatrix_jTee);

			// extracts the translation part and computes the skew symmetric matrix
			cv::Mat l_oTranslationMatrix_jtee = l_oTransformationMatrix.getTranslationMatrix()->getTranslationMat();
			cv::Mat l_oSkewSymmetricMatrix = skewSymmetricMatrix(l_oTranslationMatrix_jtee.at<double>(0, 0), l_oTranslationMatrix_jtee.at<double>(1, 0), l_oTranslationMatrix_jtee.at<double>(2, 0));

			// extracts the rotation part eeRj
			cv::Mat l_oTransformationMatrix_eeTj = l_oTransformationMatrix.getTransformationMatInverse();
			cv::Mat l_oRotationMatrix_eeRj = l_oTransformationMatrix_eeTj(cv::Rect(0, 0, 3, 3));

			// computes Tj i.e. [eeRj -eeRj*ssm; 0 eeRj]]
			cv::Mat l_oTj(6, 6, CV_64F, 0.0f);
			l_oRotationMatrix_eeRj.copyTo(l_oTj(cv::Rect(0, 0, l_oRotationMatrix_eeRj.cols, l_oRotationMatrix_eeRj.rows)));
			l_oRotationMatrix_eeRj.copyTo(l_oTj(cv::Rect(3, 3, l_oRotationMatrix_eeRj.cols, l_oRotationMatrix_eeRj.rows)));
			cv::Mat l_oRightUpperPartOfTj = -l_oRotationMatrix_eeRj * l_oSkewSymmetricMatrix;
			l_oRightUpperPartOfTj.copyTo(l_oTj(cv::Rect(3, 0, l_oRightUpperPartOfTj.cols, l_oRightUpperPartOfTj.rows)));

			// selects the right column depending on the type of joint (revolute or prismatic)
			int l_i32IndexColumn = 0;
			if (m_vJointType[l_joint] == KINEMATICS_REVOLUTE_JOINT)
				l_i32IndexColumn = 5;
			else if (m_vJointType[l_joint] == KINEMATICS_PRISMATIC_JOINT)
				l_i32IndexColumn = 2;

			// fills the numerical Jacobian
			l_oTj(cv::Rect(l_i32IndexColumn, 0, 1, l_oTj.rows)).copyTo(l_oNumericalJacobian(cv::Rect(l_i32JointIndex, 0, 1, l_oTj.rows)));

			// increments the joint index
			l_i32JointIndex++;
		}
	}

	return l_oNumericalJacobian;
}

cv::Mat Kinematics::computePseudoInverseJacobian(cv::Mat oNumericalJacobian)
{
	cv::Mat l_oPseudoInverseJacobian;
	if (false) // OpenCV SVD generates wrong pseudo inverse
	{
		// decomposes the numerical jacobian usind SVD
		cv::SVD l_oSVD;
		cv::Mat l_oW, l_oU, l_oVt;
		l_oSVD.compute(oNumericalJacobian, l_oW, l_oU, l_oVt, cv::SVD::FULL_UV);

		// converts the singular values to matrix format (OpenCV stored them in a column vector)
		cv::Mat l_oSigma(6, m_i32NbJoints, CV_64F, 0.0);
		for (int l_joint = 0; l_joint < m_i32NbJoints; l_joint++)
			l_oSigma.at<double>(l_joint, l_joint) = l_oW.at<double>(l_joint, 0);

		// computes the inverse of the sigma matrix (checks if the ratio between the current singular value and the first one is higher than a threshold)
		cv::Mat l_oSigmaInv(l_oSigma.cols, l_oSigma.rows, CV_64F, 0.0);

		for (int l_singVal = 0; l_singVal < l_oSigma.cols; l_singVal++)
		{
			if (abs(l_oSigma.at<double>(l_singVal, l_singVal) / l_oSigma.at<double>(0, 0)) < m_f64RankThreshold)
				l_oSigmaInv.at<double>(l_singVal, l_singVal) = 0.0;
			else
				l_oSigmaInv.at<double>(l_singVal, l_singVal) = 1.0 / l_oSigma.at<double>(l_singVal, l_singVal);
		}

		// computes the pseudoinverse matrix of the numerical jacobian
		l_oPseudoInverseJacobian = l_oVt * l_oSigmaInv * l_oU.t(); // = V * Sinv * Ut if U, S, V
	}
	else
	{
	
		// maps the OpenCV matrix with Eigen:
		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> oNumJacobian_Eigen(oNumericalJacobian.ptr<double>(), oNumericalJacobian.rows, oNumericalJacobian.cols);
		
		// computes the Jacobi SVD of the Numerical Jacobian
		//Eigen::JacobiSVD<Eigen::MatrixXd> svd(oNumJacobian_Eigen, Eigen::ComputeFullV | Eigen:: ComputeFullU);
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(oNumJacobian_Eigen, Eigen::ComputeFullV | Eigen:: ComputeFullU);
		svd.computeU();
		svd.computeV();
		
		// stores the U, V, S matrices
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> oUMatrix_Eigen = svd.matrixU();
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> oVMatrix_Eigen = svd.matrixV();
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> oSigmaMatrix_Eigen = svd.singularValues();
		
		// maps back to OpenCV matrices
		cv::Mat oUMatrix_OpenCV(oUMatrix_Eigen.rows(), oUMatrix_Eigen.cols(), CV_64F, oUMatrix_Eigen.data());
		cv::Mat oVMatrix_OpenCV(oVMatrix_Eigen.rows(), oVMatrix_Eigen.cols(), CV_64F, oVMatrix_Eigen.data());
		cv::Mat oSigmaMatrix_OpenCV(oSigmaMatrix_Eigen.rows(), oSigmaMatrix_Eigen.cols(), CV_64F, oSigmaMatrix_Eigen.data());

		// computes the inverse of the sigma matrix (checks if the ratio between the current singular value and the first one is higher than a threshold)
		cv::Mat l_oSigmaInv(m_i32NbJoints, 6, CV_64F, 0.0);

		for (int l_singVal = 0; l_singVal < oSigmaMatrix_OpenCV.rows; l_singVal++)
		{
			if (abs(oSigmaMatrix_OpenCV.at<double>(l_singVal, 0) / oSigmaMatrix_OpenCV.at<double>(0, 0)) < m_f64RankThreshold)
				l_oSigmaInv.at<double>(l_singVal, l_singVal) = 0.0;
			else
				l_oSigmaInv.at<double>(l_singVal, l_singVal) = 1.0 / oSigmaMatrix_OpenCV.at<double>(l_singVal, 0);
		}

		// computes the pseudo inverse of the jacobian matrix
		l_oPseudoInverseJacobian = oVMatrix_OpenCV * l_oSigmaInv * oUMatrix_OpenCV.t(); // = V * Sinv * Ut if U, S, V
	}
	
	return l_oPseudoInverseJacobian;
}


double Kinematics::computeDistanceToTarget(cv::Mat& oXtarget)
{
	// estimates the current X (position and orientation of the end-effector)
	cv::Mat l_oXcurrent_hc = computeCurrentEndEffectorPosition();
	
	// computes the Euclidian distance between the current point and the target one
	double l_fDistanceToTarget = 0.0f;
	for (int l_coord = 0; l_coord < 3; l_coord++)
	{
		l_fDistanceToTarget += pow(oXtarget.at<double>(l_coord, 0) - l_oXcurrent_hc.at<double>(l_coord, 0), 2);
	}

	l_fDistanceToTarget = sqrt(l_fDistanceToTarget);

	return l_fDistanceToTarget;
}

void Kinematics::inverseKinematics(cv::Mat oXTarget)
{
	bool l_bLoopInProgress = true;
	int l_i32NbIterations = 0;
	
	if (m_bIsDatabaseLoaded)
	{
		findClosestCandidateInDatabase(oXTarget);
	}
	else
	{
		m_vQiValues = m_vQiOriginalValues;
	}

	while (l_bLoopInProgress)
	{
		cv::Mat l_oVee = computeVee(oXTarget);
		//std::cout << "inverseKinematics : l_oVee = " << l_oVee << std::endl;
		cv::Mat l_oJNum = computeNumericalJacobian();
		//std::cout << "inverseKinematics : l_oJNum = " << l_oJNum << std::endl;
		cv::Mat l_oJPseudoInv = computePseudoInverseJacobian(l_oJNum);
		//std::cout << "inverseKinematics : l_oJPseudoInv = " << l_oJPseudoInv << std::endl;
		cv::Mat deltaQ = l_oJPseudoInv * l_oVee;
		//std::cout << "inverseKinematics : deltaQ = " << deltaQ << std::endl;
		
		int l_i32JointIndex = 0;
		for (int l_joint = 0; l_joint < m_vJointType.size(); l_joint++)
		{
			if (m_vJointType[l_joint] == KINEMATICS_REVOLUTE_JOINT || m_vJointType[l_joint] == KINEMATICS_PRISMATIC_JOINT)
			{
				m_vQiValues[l_i32JointIndex] += deltaQ.at<double>(l_i32JointIndex, 0);
				//std::cout << "(Kinematics::inverseKinematics) (joint index, joint value) = ( " << l_i32JointIndex << ", " << m_vQiValues[l_i32JointIndex]   << ")" << std::endl;
				l_i32JointIndex++;
			}
		}
		
		updateDHParameters();

		double l_fCurrentDistance = computeDistanceToTarget(oXTarget);
		std::cout << "---> #" << l_i32NbIterations << ", distance = " << l_fCurrentDistance << std::endl;
		
		l_i32NbIterations++;

		if (l_fCurrentDistance < m_f64DistanceThreshold || l_i32NbIterations > m_i32MaxNbIterations)
			l_bLoopInProgress = false;
	}
}

cv::Mat Kinematics::computeCurrentEndEffectorPosition()
{
	// estimates the current position of the end-effector)
	cv::Mat l_oTransformationMatrix_bTee = dh2FK(0);
	cv::Mat l_oXorigin_hc(4, 1, CV_64F, 0.0f); l_oXorigin_hc.at<double>(3, 0) = 1.0f; // homogeneous coordinates Xorigin = [0 0 0 1]'
	cv::Mat l_oXcurrent_hc = l_oTransformationMatrix_bTee * l_oXorigin_hc;

	return l_oXcurrent_hc;
}

void Kinematics::createKinematicDatabase(int i32NbSamples, std::string sKinematicDatabaseFilename, int i32MinValueRevoluteJoint, int i32MaxValueRevoluteJoint, int i32MinValuePrismaticJoint, int i32MaxValuePrismaticJoint)
{
	// output matrices where the random qi values and the corresponding end-effector positions will be stored
	cv::Mat l_oRandomJointValues(i32NbSamples, m_i32NbJoints, CV_64F, 0.0f);
	cv::Mat l_oRandomEndEffectorPositions(i32NbSamples, 3, CV_64F, 0.0f);

	// setups random number generators
	std::random_device l_oRandomDevice; // seed with a real random value, if available
	std::default_random_engine l_oRandomEngine(l_oRandomDevice());// create a random engine
	std::uniform_int_distribution<int> l_oUniformDistributionRevoluteJoint(i32MinValueRevoluteJoint, i32MaxValueRevoluteJoint);
	std::uniform_int_distribution<int> l_oUniformDistributionPrismaticJoint(i32MinValuePrismaticJoint, i32MaxValuePrismaticJoint);

	// loops over the samples and generates random qi values and the corresponding end-effector positions
	for (int l_sample = 0; l_sample < i32NbSamples; l_sample++)
	{
		// generates random qi values depending on the joint type i.e. revolute or prismatic
		int l_i32JointIndex = 0;
		for (int l_joint = 0; l_joint < m_vJointType.size(); l_joint++)
		{
			if (m_vJointType[l_joint] == KINEMATICS_REVOLUTE_JOINT)
			{
				int l_i32RandomRevoluteJointValue = l_oUniformDistributionRevoluteJoint(l_oRandomEngine);
				l_oRandomJointValues.at<double>(l_sample, l_i32JointIndex) = (double)l_i32RandomRevoluteJointValue;
				m_vQiValues[l_i32JointIndex] = (double)l_i32RandomRevoluteJointValue;
				l_i32JointIndex++;

			}
			else if (m_vJointType[l_joint] == KINEMATICS_PRISMATIC_JOINT)
			{
				int l_i32RandomPrismaticJointValue = l_oUniformDistributionPrismaticJoint(l_oRandomEngine);
				l_oRandomJointValues.at<double>(l_sample, l_i32JointIndex) = (double)l_i32RandomPrismaticJointValue / 100.0f;
				m_vQiValues[l_i32JointIndex] = (double)l_i32RandomPrismaticJointValue / 100.0f;
				l_i32JointIndex++;
			}
		}
		updateDHParameters();

		// computes the end-effector position
		cv::Mat l_oCurrentEndEffectorPosition = computeCurrentEndEffectorPosition();

		// stores this data to a cv::Mat
		for (int l_coord= 0; l_coord < 3; l_coord++)
			l_oRandomEndEffectorPositions.at<double>(l_sample, l_coord) = l_oCurrentEndEffectorPosition.at<double>(l_coord, 0);
		
	}
	
	// creates a file storage to store the data
	cv::FileStorage l_fsKinematicDatabaseFile(sKinematicDatabaseFilename, cv::FileStorage::WRITE);

	// checks if it was opened correctly
	if (l_fsKinematicDatabaseFile.isOpened() == false)
	{
		std::cerr << "[ERROR] (Kinematics::createKinematicDatabase) Could not open the kinematic database file for writing!";
		return;
	}

	// writes data in the file
	l_fsKinematicDatabaseFile << "randomJointValues" << l_oRandomJointValues;
	l_fsKinematicDatabaseFile << "randomEndEffectorPositions" << l_oRandomEndEffectorPositions;
	
	l_fsKinematicDatabaseFile.release();

}

void Kinematics::loadKinematicDatabase(std::string sKinematicDatabaseFilename)
{
	// checks if the filename is empty
	if (sKinematicDatabaseFilename.empty())
	{
		std::cerr << "[ERROR] (Kinematics::loadKinematicDatabase) Kinematic database filename is empty!";
		return;
	}
		
	// loads the databse file
	cv::FileStorage l_fsKinematicDatabaseFile(sKinematicDatabaseFilename, cv::FileStorage::READ);

	// checks if file was successfully opened
	if (l_fsKinematicDatabaseFile.isOpened() == false)
	{
		std::cerr << "[ERROR] (Kinematics::loadKinematicDatabase) Could not open the kinematic database file for reading!";
		return;
	}

	// reads the different variables
	cv::Mat l_oRandomJointValues;
	l_fsKinematicDatabaseFile["randomJointValues"] >> l_oRandomJointValues;
	if (l_oRandomJointValues.rows != 0 && l_oRandomJointValues.cols != 0)
		m_oRandomJointValues = l_oRandomJointValues;

	cv::Mat l_oRandomEndEffectorPositions;
	l_fsKinematicDatabaseFile["randomEndEffectorPositions"] >> l_oRandomEndEffectorPositions;
	if (l_oRandomEndEffectorPositions.rows != 0 && l_oRandomEndEffectorPositions.cols != 0)
		m_oRandomEndEffectorPositions = l_oRandomEndEffectorPositions;

	

	m_bIsDatabaseLoaded = true;
}

void Kinematics::loadDHParameters(std::string sDHParametersFilename)
{
	// checks if the filename is empty
	if (sDHParametersFilename.empty())
	{
		std::cerr << "[ERROR] (Kinematics::loadDHParameters) DH Parameters filename is empty!";
		return;
	}
		
	// loads the databse file
	cv::FileStorage l_fsDHParametersFile(sDHParametersFilename, cv::FileStorage::READ);
		
	// checks if file was successfully opened
	if (l_fsDHParametersFile.isOpened() == false)
	{
		std::cerr << "[ERROR] (Kinematics::loadDHParameters) Could not open the DH Parameters file for reading!";
		return;
	}

	// reads the different variables
	std::vector<double> l_vDHTheta;
	l_fsDHParametersFile["dhThetaParams"] >> l_vDHTheta;
	if (l_vDHTheta.size() > 0)
		m_vDHTheta = l_vDHTheta;
	
	std::vector<double> l_vDHD;
	l_fsDHParametersFile["dhDParams"] >> l_vDHD;
	if (l_vDHD.size() > 0)
		m_vDHD = l_vDHD;
	
	std::vector<double> l_vDHA;
	l_fsDHParametersFile["dhAParams"] >> l_vDHA;
	if (l_vDHA.size() > 0)
		m_vDHA = l_vDHA;
	
	std::vector<double> l_vDHAlpha;
	l_fsDHParametersFile["dhAlphaParams"] >> l_vDHAlpha;
	if (l_vDHAlpha.size() > 0)
		m_vDHAlpha = l_vDHAlpha;
	
	std::vector<int> l_vJointType;
	l_fsDHParametersFile["jointType"] >> l_vJointType;
	if (l_vJointType.size() > 0)
		m_vJointType = l_vJointType;
	
	l_fsDHParametersFile["nbJoints"] >> m_i32NbJoints;

	std::vector<double> l_vJointOffset;
	l_fsDHParametersFile["jointOffset"] >> l_vJointOffset;
	if (l_vJointOffset.size() > 0)
		m_vJointOffset = l_vJointOffset;
	
	std::vector<double> l_vQiValues;
	l_fsDHParametersFile["qiValues"] >> l_vQiValues;
	if (l_vQiValues.size() > 0)
	{
		m_vQiValues = l_vQiValues;
		m_vQiOriginalValues = l_vQiValues;
	}
		
	updateDHParameters();
	
	m_bAreDHParametersSet = true;
}


void Kinematics::loadIKParameters(std::string sIKParametersFilename)
{
	// checks if the filename is empty
	if (sIKParametersFilename.empty())
	{
		std::cerr << "[ERROR] (Kinematics::loadIKParameters) IK Parameters filename is empty!";
		return;
	}
		
	// loads the databse file
	cv::FileStorage l_fsIKParametersFile(sIKParametersFilename, cv::FileStorage::READ);

	// checks if file was successfully opened
	if (l_fsIKParametersFile.isOpened() == false)
	{
		std::cerr << "[ERROR] (Kinematics::loadIKParameters) Could not open the IK Parameters file for reading!";
		return;
	}

	// reads the different variables
	l_fsIKParametersFile["maxNbIterations"] >> m_i32MaxNbIterations;
	l_fsIKParametersFile["pController"] >> m_f64PController;
	l_fsIKParametersFile["rankThreshold"] >> m_f64RankThreshold;
	l_fsIKParametersFile["distanceThreshold"] >> m_f64DistanceThreshold;

	m_bAreDHParametersSet = true;
}


void Kinematics::findClosestCandidateInDatabase(cv::Mat& oXtarget)
{
	double l_fClosestSampleDistanceToTarget = std::numeric_limits<double>::max();
	int l_i32ClosestSampleIndex = 0;

	//std::cout << "m_oRandomEndEffectorPositions.rows = " << m_oRandomEndEffectorPositions.rows << std::endl;
	//std::cout << "m_oRandomEndEffectorPositions.cols = " << m_oRandomEndEffectorPositions.cols << std::endl;
	// looks for the closest candidate in the database
	for (int l_sample = 0; l_sample < m_oRandomEndEffectorPositions.rows; l_sample++)
	{
		double l_fCurrentSampleDistanceToTarget = 0.0f;
		for (int l_coord = 0; l_coord < 3; l_coord++)
			l_fCurrentSampleDistanceToTarget += pow( oXtarget.at<double>(l_coord, 0) - m_oRandomEndEffectorPositions.at<double>(l_sample, l_coord), 2);

		if (l_fCurrentSampleDistanceToTarget < l_fClosestSampleDistanceToTarget)
		{
			l_fClosestSampleDistanceToTarget = l_fCurrentSampleDistanceToTarget;
			l_i32ClosestSampleIndex = l_sample;
		}
	}
	
	std::cout << "l_fClosestSampleDistanceToTarget = " << l_fClosestSampleDistanceToTarget << std::endl;
	std::cout << "l_i32ClosestSampleIndex = " << l_i32ClosestSampleIndex << std::endl;

	// change the DH parameters to become the closest candidate in the database
	int l_i32JointIndex = 0;
	for (int l_joint = 0; l_joint < m_vJointType.size(); l_joint++)
	{
		if (m_vJointType[l_joint] == KINEMATICS_REVOLUTE_JOINT || m_vJointType[l_joint] == KINEMATICS_PRISMATIC_JOINT)
		{
			m_vQiValues[l_joint] = m_oRandomJointValues.at<double>(l_i32ClosestSampleIndex, l_i32JointIndex);
			l_i32JointIndex++;
		}
	}
	
	updateDHParameters();
}
