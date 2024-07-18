#include "TranslationMatrix.h"


TranslationMatrix::TranslationMatrix():
m_bIsMatrixComputed(false), 
m_fTx(0.0f), m_fTy(0.0f), m_fTz(0.0f)
{}

TranslationMatrix::TranslationMatrix(double fTx, double fTy, double fTz):
m_bIsMatrixComputed(false), 
m_fTx(fTx), m_fTy(fTy), m_fTz(fTz)
{
	compute();
}

TranslationMatrix::TranslationMatrix(cv::Mat& oTranslationMatrix)
{
	m_oTranslationMatrix = oTranslationMatrix;

	m_fTx = oTranslationMatrix.at<double>(0, 0);
	m_fTy = oTranslationMatrix.at<double>(1, 0);
	m_fTz = oTranslationMatrix.at<double>(2, 0);

	m_bIsMatrixComputed = true;
}

TranslationMatrix::~TranslationMatrix()
{
}

void TranslationMatrix::setTranslationComponents(double fTx, double fTy, double fTz)
{
	m_fTx = fTx;
	m_fTy = fTy;
	m_fTz = fTz;
	
	m_bIsMatrixComputed = false;
}

cv::Mat TranslationMatrix::getTranslationMat()
{
	return m_oTranslationMatrix;
}

void TranslationMatrix::compute()
{
	// creates matrices to store individual roation matrices Rx, Ry and Rz before product
	cv::Mat l_oTranslationMatrix(3, 1, CV_64F);
	l_oTranslationMatrix.at<double>(0, 0) = m_fTx; 	
	l_oTranslationMatrix.at<double>(1, 0) = m_fTy; 	
	l_oTranslationMatrix.at<double>(2, 0) = m_fTz; 	
	
	m_oTranslationMatrix = l_oTranslationMatrix;
	
	m_bIsMatrixComputed = true;

}

void TranslationMatrix::disp()
{
	if (m_bIsMatrixComputed)
	{
		std::cout << "Translation matrix =" << std::endl;
		std::cout << m_oTranslationMatrix << std::endl;
	}
	else
		std::cerr << "(TranslationMatrix::disp) The translation matrix has not been  computed yet!" << std::endl;
}


