#ifndef TRANSFORMATION_MATRIX_
#define TRANSFORMATION_MATRIX_

#define _USE_MATH_DEFINES
#include <cmath>


#include "opencv2/opencv.hpp"

#include "RotationMatrix.h"
#include "TranslationMatrix.h"

class TransformationMatrix
{
	public:
		TransformationMatrix();
	
		TransformationMatrix(double fRoll, double fPitch, double fYaw, int i32Order, double fTx, double fTy, double fTz);

		TransformationMatrix(cv::Mat& oTransformationMatrix);
	
		TransformationMatrix(cv::Mat& rotationMatrix, cv::Mat& translationMatrix);
	
		~TransformationMatrix();
	
		void setRotationComponents(double fRoll, double fPitch, double fYaw, int i32Order);
	
		void setTranslationComponents(double fTx, double fTy, double fTz);

		cv::Mat getTransformationMat();

		cv::Mat getTransformationMatInverse();

		RotationMatrix* getRotationMatrix();

		TranslationMatrix* getTranslationMatrix();

		/**
		* \brief Computes the transformation matrix from the roll, pitch and yaw angles, the order of operations and the translation components
		*/
		void compute();
		
		/**
		* \brief Computes the inverse of the transformation matrix i.e. 	(R^-1	|-R^-1t)
		*											-----------------
		*											(   0		|   1)
		*/
		void inverse();
		
		/**
		* \brief Displays the transformation matrix and its inverse (if computed)
		*/
		void disp();

	
	private:
		bool m_bAreRotationComponentsSet;
		bool m_bAreTranslationComponentsSet;
		bool m_bIsMatrixComputed;		/**< boolean value to check if the rotation matrix is already computed */
		bool m_bIsMatrixInverseComputed;	/**< boolean value to check if the rotation matrix inverse is already computed */
		RotationMatrix* m_pRotationMatrix;	/**< pointer to a rotation matrix (3x3) */
		TranslationMatrix* m_pTranslationMatrix;	/**< pointer to a translation matrix (3x1) */
		cv::Mat m_oTransformationMatrix;
		cv::Mat m_oTransformationMatrixInverse;
		
};


#endif