#ifndef ROTATION_MATRIX_
#define ROTATION_MATRIX_

#define _USE_MATH_DEFINES
#include <cmath>

#include "opencv2/opencv.hpp"

#define RxRyRz 0
#define RzRyRx 1

class RotationMatrix
{
	public:
		RotationMatrix();
	
		RotationMatrix(double fRoll, double fPitch, double fYaw, int i32Order);
	
		RotationMatrix(cv::Mat& oRotationMatrix);

		~RotationMatrix();

		/**
		* \brief Sets new roll, pitch and yaw angles.
		* \param [in] roll : value of the roll angles in degrees
		* \param [in] pitch : value of the pitch angles in degrees
		* \param [in] yaw : value of the yaw angles in degrees
		*/
		void setRollPitchYawAngles(double fRoll, double fPitch, double fYaw);
	
		/**
		* \brief Sets the order of operation RxRyRz or RzRyRx.
		* \param [in] order : order of operation (either RxRyRz or RzRyRx i.e. Cardanian representation) 
		*/
		void setOrderOfOperations(int i32Order);
	
		cv::Mat getRotationMat();
		
		cv::Mat getRotationMatInverse();
		
		/**
		* \brief Computes the rotation matrix from the roll, pitch and yaw angles and the order of operations
		*/
		void compute();
		
		/**
		* \brief Computes the inverse of the rotation matrix i.e. its transpose
		*/
		void inverse();
		
		/**
		* \brief Displays the rotation matrix and its inverse (if computed)
		*/
		void disp();

	private:
		/**
		* \brief Converts an angle from degrees to radians
		* \param [in] angleInDegree : input angle in degrees to be converted in radians
		*\return a double corresponding to the angle in radians
		*/
		double deg2rad(double angleInDegrees);
	
		/**
		* \brief Converts an angle from radian to degree
		* \param [in] angleInRadians : input angle in radians to be converted in degrees
		*\return a double corresponding to the angle in degrees
		*/
		double rad2deg(double angleInRadians);

	private:
		bool m_bIsMatrixComputed;		/**< boolean value to check if the rotation matrix is already computed */
		bool m_bIsMatrixInverseComputed;	/**< boolean value to check if the rotation matrix inverse is already computed */
		double m_fRollAngle;				/**< roll angle in degrees */
		double m_fPitchAngle;				/**< pitch angle in degrees */
		double m_fYawAngle;				/**< yaw angle in degrees */
		int m_i32OrderOfOperations;		/**< order of operations */
		cv::Mat m_oRotationMatrix; 		/**< rotation matrix */
		cv::Mat m_oRotationMatrixInverse; 	/**< inverse rotation matrix */

};


#endif