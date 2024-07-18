#ifndef KINEMATICS_
#define KINEMATICS_

#include <vector>
#include <string>
#include <random>
#include <limits>

#include "opencv2/opencv.hpp"

#include "Dense"

#include "TransformationMatrix.h"

#define KINEMATICS_REVOLUTE_JOINT 1
#define KINEMATICS_PRISMATIC_JOINT 2

class Kinematics
{
	public:
		// Constructors
		Kinematics();
		Kinematics(std::vector<double>& vDHTheta, std::vector<double>& vDHD, std::vector<double>& vDHA, std::vector<double>& vDHAlpha, std::vector<int>& vJointType, std::vector<double>& vJointOffset, std::vector<double>& vQiValues, int i32NbJoints, int i32MaxNbIterations, double f64PController, double f64RankThreshold, double f64DistanceThreshold);
		// Destructor
		~Kinematics();
			
		// DH Params
		void loadDHParameters(std::string sDHParametersFilename);
		void setDHParameters(std::vector<double>& vDHTheta, std::vector<double>& vDHD, std::vector<double>& vDHA, std::vector<double>& vDHAlpha, std::vector<int>& vJointType, std::vector<double>& vJointOffset, std::vector<double>& vQiValues, int i32NbJoints);
		void updateDHParameters(std::vector<double> vQiValues);
		void updateDHParameters();
		std::vector<double> getQiValues();

		// FK
		cv::Mat dh2FK(int i32StartJointNumber);
		cv::Mat computeCurrentEndEffectorPosition();
	
		// IK
		void loadIKParameters(std::string sIKParametersFilename);
		void inverseKinematics(cv::Mat oXTarget);

		// Database
		void createKinematicDatabase(int i32NbSamples, std::string sKinematicDatabaseFilename, int i32MinValueRevoluteJoint, int i32MaxValueRevoluteJoint, int i32MinValuePrismaticJoint, int i32MaxValuePrismaticJoint);
		void loadKinematicDatabase(std::string sKinematicDatabaseFilename);
		void findClosestCandidateInDatabase(cv::Mat& oXtarget);

	private:
		cv::Mat skewSymmetricMatrix(double fX, double fY, double fZ);

		cv::Mat computeVee(cv::Mat oXtarget);

		cv::Mat computeNumericalJacobian();

		cv::Mat computePseudoInverseJacobian(cv::Mat oNumericalJacobian);

		double computeDistanceToTarget(cv::Mat& oXtarget);

	private:
		// DH Params
		bool m_bAreDHParametersSet;
		int m_i32NbJoints;
		std::vector<double> m_vDHTheta;
		std::vector<double> m_vDHD;
		std::vector<double> m_vDHA;
		std::vector<double> m_vDHAlpha;
		std::vector<double> m_vJointOffset;
		std::vector<int> m_vJointType;
		std::vector<double> m_vQiValues;
		std::vector<double> m_vQiOriginalValues;
		// Database
		bool m_bIsDatabaseLoaded;
		cv::Mat m_oRandomEndEffectorPositions;
		cv::Mat m_oRandomJointValues;
		// IK
		int m_i32MaxNbIterations;
		double m_f64PController;
		double m_f64RankThreshold;
		double m_f64DistanceThreshold;
};

#endif