#include "TransformationMatrix.h"


TransformationMatrix::TransformationMatrix():
m_bAreRotationComponentsSet(false), m_bAreTranslationComponentsSet(false),
m_bIsMatrixComputed(false), m_bIsMatrixInverseComputed(false)
{
	m_pRotationMatrix = new RotationMatrix();
	m_pTranslationMatrix = new TranslationMatrix();
}

TransformationMatrix::TransformationMatrix(double fRoll, double fPitch, double fYaw, int i32Order, double fTx, double fTy, double fTz):
m_bAreRotationComponentsSet(true), m_bAreTranslationComponentsSet(true),
m_bIsMatrixComputed(false), m_bIsMatrixInverseComputed(false)
{
	m_pRotationMatrix = new RotationMatrix(fRoll, fPitch, fYaw, i32Order);
	m_pTranslationMatrix = new TranslationMatrix(fTx, fTy, fTz);
	
	compute();
}

TransformationMatrix::TransformationMatrix(cv::Mat& oRotationMatrix, cv::Mat& oTranslationMatrix) :
	m_bAreRotationComponentsSet(false), m_bAreTranslationComponentsSet(false),
	m_bIsMatrixComputed(true), m_bIsMatrixInverseComputed(false)
{
	// feeds the translation part
	m_pTranslationMatrix = new TranslationMatrix(oTranslationMatrix);

	// feeds the rotation part
	m_pRotationMatrix = new RotationMatrix(oRotationMatrix);

	// feeds the transformation matrix
	cv::Mat l_oTransformationMatrix(4, 4, CV_64F, 0.0f);
	oTranslationMatrix.copyTo(l_oTransformationMatrix(cv::Rect(3, 0, oTranslationMatrix.cols, oTranslationMatrix.rows)));
	oRotationMatrix.copyTo(l_oTransformationMatrix(cv::Rect(0, 0, oRotationMatrix.cols, oRotationMatrix.rows)));
	l_oTransformationMatrix.at<double>(3, 3) = 1.0f;
	m_oTransformationMatrix = l_oTransformationMatrix;

	inverse();
}

TransformationMatrix::TransformationMatrix(cv::Mat& oTransformationMatrix) :
	m_bAreRotationComponentsSet(false), m_bAreTranslationComponentsSet(false),
	m_bIsMatrixComputed(true), m_bIsMatrixInverseComputed(false),
	m_oTransformationMatrix(oTransformationMatrix)
{
	// feeds the translation part
	cv::Mat l_oTranslationMatrix(3, 1, CV_64F);
	oTransformationMatrix(cv::Rect(3, 0, 1, 3)).copyTo(l_oTranslationMatrix);
	m_pTranslationMatrix = new TranslationMatrix(l_oTranslationMatrix);
	
	// feeds the rotation part
	cv::Mat l_oRotationMatrix(3, 3, CV_64F);
	oTransformationMatrix(cv::Rect(0, 0, 3, 3)).copyTo(l_oRotationMatrix);
	m_pRotationMatrix = new RotationMatrix(l_oRotationMatrix);
		
	inverse();
}

TransformationMatrix::~TransformationMatrix()
{
	// lambda function to delete and nullify any pointer
	auto deleteAndNullify = [](auto pointer) -> void
	{
		if (nullptr != pointer)
		{
			delete pointer;
			pointer = nullptr;
		}
	};
	
	// release allocated memory
	deleteAndNullify(m_pRotationMatrix);
	deleteAndNullify(m_pTranslationMatrix);
}

void TransformationMatrix::setRotationComponents(double fRoll, double fPitch, double fYaw, int i32Order)
{
	if (nullptr != m_pRotationMatrix)
	{
		delete m_pRotationMatrix;
		m_pRotationMatrix = nullptr;
	}
	
	m_pRotationMatrix = new RotationMatrix(fRoll, fPitch, fYaw, i32Order);
	
	m_bAreRotationComponentsSet = true;
	m_bIsMatrixComputed = false;
	m_bIsMatrixInverseComputed = false;
}


void TransformationMatrix::setTranslationComponents(double fTx, double fTy, double fTz)
{
	if (nullptr != m_pTranslationMatrix)
	{
		delete m_pTranslationMatrix;
		m_pTranslationMatrix = nullptr;
	}
	
	m_pTranslationMatrix = new TranslationMatrix(fTx, fTy, fTz);
	
	m_bAreTranslationComponentsSet = true;
	m_bIsMatrixComputed = false;
	m_bIsMatrixInverseComputed = false;
}


cv::Mat TransformationMatrix::getTransformationMat()
{
	if (!m_bIsMatrixComputed)
		std::cerr << "(TransformationMatrix::getTransformationMatrix) Transformation matrix has not been computed yet!" << std::endl;

	return m_oTransformationMatrix;
}

cv::Mat TransformationMatrix::getTransformationMatInverse()
{
	if (!m_bIsMatrixInverseComputed)
		std::cerr << "(TransformationMatrix::getTransformationMatrixInverse) Transformation matrix inverse has not been computed yet!" << std::endl;

	return m_oTransformationMatrixInverse;
}

RotationMatrix* TransformationMatrix::getRotationMatrix()
{
	if (!m_bIsMatrixComputed)
		std::cerr << "(TransformationMatrix::getRotationMatrix) Transformation matrix has not been computed yet!" << std::endl;

	return m_pRotationMatrix;
}

TranslationMatrix* TransformationMatrix::getTranslationMatrix()
{
	if (!m_bIsMatrixComputed)
		std::cerr << "(TransformationMatrix::getTranslationMatrix) Transformation matrix has not been computed yet!" << std::endl;

	return m_pTranslationMatrix;
}

void TransformationMatrix::compute()
{
	if (m_bAreRotationComponentsSet && m_bAreTranslationComponentsSet)
	{
		cv::Mat l_oTransformationMatrix(4, 4, CV_64F); // homogeneous coordinates
		/*	(R (3x3)	|	t (3x1) 	)
		*	-----------------------------------
		*	(   0	(1x3)	|	1 (1x1)	)
		*/
		
		// retrieves the rotation components and copies them into the transformation matrix
		cv::Mat l_oRotationMatrix = m_pRotationMatrix->getRotationMat();
		l_oRotationMatrix.copyTo(l_oTransformationMatrix(cv::Rect(0, 0, l_oRotationMatrix.cols, l_oRotationMatrix.rows)));
		
		// retrieves the translation components and copies them into the transformation matrix
		cv::Mat l_oTranslationMatrix = m_pTranslationMatrix->getTranslationMat();
		l_oTranslationMatrix.copyTo(l_oTransformationMatrix(cv::Rect(3, 0, l_oTranslationMatrix.cols, l_oTranslationMatrix.rows)));
		
		// adds the homogeneous coordinate part
		l_oTransformationMatrix.at<double>(3, 0) = 0.0f;
		l_oTransformationMatrix.at<double>(3, 1) = 0.0f;
		l_oTransformationMatrix.at<double>(3, 2) = 0.0f;
		l_oTransformationMatrix.at<double>(3, 3) = 1.0f;
		
		// stores it in the member variable
		m_oTransformationMatrix = l_oTransformationMatrix;
		
		m_bIsMatrixComputed = true;
	}
	else
	{
		std::cerr << "(TransformationMatrix::compute) Rotation and/or translation components are not set!" << std::endl;
	}

}

void TransformationMatrix::inverse()
{
	if (m_bIsMatrixComputed)
	{
		cv::Mat l_oTransformationMatrixInverse(4, 4, CV_64F); // homogeneous coordinates
		/*	(R^-1 (3x3)	|	-R^-1t (3x1) )
		*	-----------------------------------
		*	(   0	(1x3)		|	1 (1x1)	)
		*/
		
		// retrieves the rotation components and copies them into the transformation matrix inverse
		cv::Mat l_oRotationMatrixInverse = m_pRotationMatrix->getRotationMatInverse();
		l_oRotationMatrixInverse.copyTo(l_oTransformationMatrixInverse(cv::Rect(0, 0, l_oRotationMatrixInverse.cols, l_oRotationMatrixInverse.rows)));
		
		// retrieves the translation components 
		cv::Mat l_oTranslationMatrix = m_pTranslationMatrix->getTranslationMat();
		// computes the inverse translation components 
		cv::Mat l_oTranslationMatrixInverse = -l_oRotationMatrixInverse * l_oTranslationMatrix;
		// copies them into the transformation matrix inverse
		l_oTranslationMatrixInverse.copyTo(l_oTransformationMatrixInverse(cv::Rect(3, 0, l_oTranslationMatrixInverse.cols, l_oTranslationMatrixInverse.rows)));
		
		// adds the homogeneous coordinate part
		l_oTransformationMatrixInverse.at<double>(3, 0) = 0.0f;
		l_oTransformationMatrixInverse.at<double>(3, 1) = 0.0f;
		l_oTransformationMatrixInverse.at<double>(3, 2) = 0.0f;
		l_oTransformationMatrixInverse.at<double>(3, 3) = 1.0f;
		
		// stores it in the member variable
		m_oTransformationMatrixInverse = l_oTransformationMatrixInverse;
		m_bIsMatrixInverseComputed =true;
	}
	else
		std::cerr << "(TransformationMatrix::inverse) The rotation matrix has not been  computed yet!" << std::endl;
}

void TransformationMatrix::disp()
{
	if (m_bIsMatrixComputed)
	{
		std::cout << "Transformation matrix =" << std::endl;
		std::cout << m_oTransformationMatrix << std::endl;
	}
	else
		std::cerr << "(TransformationMatrix::disp) The rotation matrix has not been  computed yet!" << std::endl;
	
	if (m_bIsMatrixInverseComputed)
	{	
		std::cout << "Transformation matrix Inverse =" << std::endl;
		std::cout << m_oTransformationMatrixInverse << std::endl;
	}
	else
		std::cerr << "(TransformationMatrix::disp) The rotation matrix inverse has not been  computed yet!" << std::endl;
	
}


