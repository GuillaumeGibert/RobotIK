#include "RotationMatrix.h"


RotationMatrix::RotationMatrix():
m_bIsMatrixComputed(false), m_bIsMatrixInverseComputed(false),
m_fRollAngle(0.0f), m_fPitchAngle(0.0f), m_fYawAngle(0.0f),
m_i32OrderOfOperations(-1)
{}

RotationMatrix::RotationMatrix(double fRoll, double fPitch, double fYaw, int i32Order):
m_bIsMatrixComputed(false), m_bIsMatrixInverseComputed(false),
m_fRollAngle(fRoll), m_fPitchAngle(fPitch), m_fYawAngle(fYaw),
m_i32OrderOfOperations(i32Order)
{
	compute();
	inverse();
}

RotationMatrix::RotationMatrix(cv::Mat& oRotationMatrix) :
	m_bIsMatrixComputed(true), m_bIsMatrixInverseComputed(false),
	m_fRollAngle(0.0f), m_fPitchAngle(0.0f), m_fYawAngle(0.0f),
	m_i32OrderOfOperations(-1),
	m_oRotationMatrix(oRotationMatrix)
{
	inverse();
}


RotationMatrix::~RotationMatrix()
{
}

void RotationMatrix::setRollPitchYawAngles(double fRoll, double fPitch, double fYaw)
{
	m_fRollAngle = fRoll;
	m_fPitchAngle = fPitch;
	m_fYawAngle = fYaw;
	
	m_bIsMatrixComputed = false;
}

void RotationMatrix::setOrderOfOperations(int i32Order)
{
	m_i32OrderOfOperations = i32Order;
	
	m_bIsMatrixComputed = false;
}

cv::Mat RotationMatrix::getRotationMat()
{
	if (!m_bIsMatrixComputed)
		std::cerr << "(RotationMatrix::getRotationMatrix) Rotation matrix has not been computed yet!" << std::endl;

	return m_oRotationMatrix;
}

cv::Mat RotationMatrix::getRotationMatInverse()
{
	if (!m_bIsMatrixInverseComputed)
		std::cerr << "(RotationMatrix::getRotationMatrix) Rotation matrix inverse has not been computed yet!" << std::endl;

	return m_oRotationMatrixInverse;
}

double RotationMatrix::deg2rad(double angleInDegrees)
{
	return (double)(angleInDegrees*M_PI/180.0f);
}

double RotationMatrix::rad2deg(double angleInRadians)
{
	return (double)(angleInRadians/M_PI*180.0f);
}

void RotationMatrix::compute()
{
	double thetaX, thetaY, thetaZ = 0.0f;
	
	// computes the theta angles depending on the order of operations
	switch(m_i32OrderOfOperations)
	{
		case RxRyRz:
			thetaX = deg2rad(m_fYawAngle);
			thetaY = deg2rad(m_fPitchAngle);
			thetaZ = deg2rad(m_fRollAngle);
			break;
		case RzRyRx:
			thetaX = deg2rad(m_fRollAngle);
			thetaY = deg2rad(m_fPitchAngle);
			thetaZ = deg2rad(m_fYawAngle);
			break;
		default:
			std::cerr << "(RotationMatrix::compute) The order of operation is unknown!" << std::endl;
			return;
	}
	
	// creates matrices to store individual roation matrices Rx, Ry and Rz before product
	cv::Mat l_oRx(3, 3, CV_64F);
	l_oRx.at<double>(0, 0) = 1; 	l_oRx.at<double>(0, 1) = 0; 			l_oRx.at<double>(0, 2) = 0;
	l_oRx.at<double>(1, 0) = 0; 	l_oRx.at<double>(1, 1) = cos(thetaX);	l_oRx.at<double>(1, 2) = -sin(thetaX);
	l_oRx.at<double>(2, 0) = 0; 	l_oRx.at<double>(2, 1) = sin(thetaX); 	l_oRx.at<double>(2, 2) = cos(thetaX);
	
	cv::Mat l_oRy(3, 3, CV_64F);
	l_oRy.at<double>(0, 0) = cos(thetaY); 	l_oRy.at<double>(0, 1) = 0;		l_oRy.at<double>(0, 2) = sin(thetaY);
	l_oRy.at<double>(1, 0) = 0; 			l_oRy.at<double>(1, 1) = 1; 		l_oRy.at<double>(1, 2) = 0;
	l_oRy.at<double>(2, 0) = -sin(thetaY); 	l_oRy.at<double>(2, 1) = 0; 		l_oRy.at<double>(2, 2) = cos(thetaY);
	
	cv::Mat l_oRz(3, 3, CV_64F);
	l_oRz.at<double>(0, 0) = cos(thetaZ); 	l_oRz.at<double>(0, 1) = -sin(thetaZ); 	l_oRz.at<double>(0, 2) = 0;
	l_oRz.at<double>(1, 0) = sin(thetaZ);		l_oRz.at<double>(1, 1) = cos(thetaZ); 	l_oRz.at<double>(1, 2) = 0;
	l_oRz.at<double>(2, 0) = 0; 			l_oRz.at<double>(2, 1) = 0; 			l_oRz.at<double>(2, 2) = 1;
	
	// computes the rotation matrix
	switch(m_i32OrderOfOperations)
	{
		case RxRyRz:
			m_oRotationMatrix = l_oRx * l_oRy * l_oRz;
			break;
		case RzRyRx:
			m_oRotationMatrix = l_oRz * l_oRy * l_oRx;
			break;
		default:
			std::cerr << "(RotationMatrix::compute) The order of operation is unknown!" << std::endl;
			return;
	}
	
	m_bIsMatrixComputed = true;

}

void RotationMatrix::inverse()
{
	if (m_bIsMatrixComputed)
	{
		m_oRotationMatrixInverse = m_oRotationMatrix.t();
		m_bIsMatrixInverseComputed =true;
	}
	else
		std::cerr << "(RotationMatrix::inverse) The rotation matrix has not been  computed yet!" << std::endl;
}

void RotationMatrix::disp()
{
	if (m_bIsMatrixComputed)
	{
		std::cout << "Rotation matrix =" << std::endl;
		std::cout << m_oRotationMatrix << std::endl;
	}
	else
		std::cerr << "(RotationMatrix::disp) The rotation matrix has not been  computed yet!" << std::endl;
	
	if (m_bIsMatrixInverseComputed)
	{	
		std::cout << "Rotation matrix Inverse =" << std::endl;
		std::cout << m_oRotationMatrixInverse << std::endl;
	}
	else
		std::cerr << "(RotationMatrix::disp) The rotation matrix inverse has not been  computed yet!" << std::endl;
	
}


