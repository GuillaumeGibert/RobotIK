#ifndef TRANSLATION_MATRIX_
#define TRANSLATION_MATRIX_

#include "opencv2/opencv.hpp"


class TranslationMatrix
{
	public:
		TranslationMatrix();
	
		TranslationMatrix(double fTx, double fTy, double fTz);

		TranslationMatrix(cv::Mat& oTranslationMatrix);
	
		~TranslationMatrix();

		/**
		* \brief Sets new translation components.
		* \param [in] fTx : value of the translation along x-axis in meters
		* \param [in] fTy : value of the translation along y-axis in meters
		* \param [in] fTz : value of the translation along z-axis in meters
		*/
		void setTranslationComponents(double fTx, double fTy, double fTz);
	
		cv::Mat getTranslationMat();
	
		/**
		* \brief Computes the rotation matrix from the roll, pitch and yaw angles and the order of operations
		*/
		void compute();
		
		/**
		* \brief Displays the translation matrix(if computed)
		*/
		void disp();

	private:
		bool m_bIsMatrixComputed;		/**< boolean value to check if the translation matrix is already computed */
		double m_fTx;					/**< translation along x in meters */
		double m_fTy;					/**< translation along y in meters */
		double m_fTz;					/**< translation along z in meters */
		cv::Mat m_oTranslationMatrix; 		/**< translation matrix */
		
};


#endif