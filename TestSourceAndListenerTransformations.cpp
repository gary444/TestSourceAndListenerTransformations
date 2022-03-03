// TestSourceAndListenerTransformations.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

typedef Eigen::Matrix<float, 4, 4, Eigen::RowMajor> Mat4f_RM;
typedef Eigen::Matrix<float, 4, 4, Eigen::ColMajor> Mat4f_CM;
typedef Eigen::Matrix<float, 3, 3, Eigen::ColMajor> Mat3f_CM;


struct SpatialData
{
	float lst_to_src_azimuth = 0.f;
	float lst_to_src_elevation = 0.f;
	float src_to_lst_azimuth = 0.f;
	float src_to_lst_elevation = 0.f;
	float lst_to_src_dist = 0.f;
	Mat4f_CM lst_transform;
	Mat4f_CM src_transform;

	float ear_L_to_mouth_azimuth = 0.f;
	float ear_L_to_mouth_elevation = 0.f;
	float ear_R_to_mouth_azimuth = 0.f;
	float ear_R_to_mouth_elevation = 0.f;
	float mouth_to_ear_L_azimuth = 0.f;
	float mouth_to_ear_L_elevation = 0.f;
	float mouth_to_ear_R_azimuth = 0.f;
	float mouth_to_ear_R_elevation = 0.f;

	float ear_L_to_mouth_dist = 0.f;
	float ear_R_to_mouth_dist = 0.f;

};

// Ear and mouth offsets from head transform
const Mat4f_CM  mouth_offset = Eigen::Affine3f(Eigen::Translation3f(0.f, -0.11f, -0.09f)).matrix();
const Mat4f_CM  ear_offset_L = Eigen::Affine3f(Eigen::Translation3f(0.11f, -0.028f, 0.02f)).matrix();
const Mat4f_CM  ear_offset_R = Eigen::Affine3f(Eigen::Translation3f(-0.11f, -0.028f, 0.02f)).matrix();

inline float FastClip(float x, float minval, float maxval) { return (fabsf(x - minval) - fabsf(x - maxval) + (minval + maxval)) * 0.5f; }


void getAzimuthAndElevationFromVector(const Eigen::Vector4f& vec, float& azimuth, float& elevation) {

	azimuth = (fabsf(vec[2]) < 0.001f) ? 0.0f : atan2f(vec[0], vec[2]);
	if (azimuth < 0.0f)
		azimuth += 2.0f * M_PI;
	azimuth = FastClip(azimuth * (180.f / M_PI), 0.0f, 360.0f);

	elevation = atan2f(vec[1], sqrtf(vec[0] * vec[0] + vec[2] * vec[2]) + 0.001f) * (180.f / M_PI);
}
void getAzimuthAndElevationFromVector(float dir_x, float dir_y, float dir_z, float& azimuth, float& elevation) {

	Eigen::Vector4f vec{ dir_x, dir_y, dir_z, 1.f };
	getAzimuthAndElevationFromVector(vec, azimuth, elevation);
}

int main()
{
	using namespace Eigen;

	SpatialData spatialData;

	Vector4f p = { 0,0,0,1 };

	const float a = 90.f * (M_PI / 180.f); // rotation angle
	const Vector3f  translate = { 0,1.5,-2 };

	Mat4f_CM lst, src;
	{
		Transform<float, 3, Eigen::Affine> t;
		t = Translation<float, 3>(translate);
		t.rotate(AngleAxis<float>(a, Vector3f::UnitY()));
		lst = t.matrix();
	}
	{
		Transform<float, 3, Eigen::Affine> t;
		t = Translation<float, 3>(Vector3f(0,0,0));
		t.rotate(AngleAxis<float>(M_PI, Vector3f::UnitY()));
		src = t.matrix();
	}

	Vector4f lst_wp = lst.col(3);
	Vector4f src_wp = src.col(3);

	std::cout << "Listener matrix = \n" << lst << std::endl;
	std::cout << "Listener position = \n" << lst_wp << std::endl;

	const Mat4f_CM ear_transform_L = lst * ear_offset_L;
	const Mat4f_CM ear_transform_R = lst * ear_offset_R;
	const Mat4f_CM mouth_transform = src * mouth_offset;

	const Eigen::Vector4f ear_L_wp = ear_transform_L.col(3);
	const Eigen::Vector4f ear_R_wp = ear_transform_R.col(3);
	const Eigen::Vector4f mouth_wp = mouth_transform.col(3);

	std::cout << "ear_L_wp = \n" << ear_L_wp << std::endl;
	std::cout << "ear_R_wp = \n" << ear_R_wp << std::endl;
	std::cout << "mouth_wp = \n" << mouth_wp << std::endl;

	const Eigen::Vector4f mouth_in_ear_L_coords = ear_transform_L.inverse() * mouth_wp;
	const Eigen::Vector4f mouth_in_ear_R_coords = ear_transform_R.inverse() * mouth_wp;
	const Eigen::Vector4f ear_L_in_mouth_coords = mouth_transform.inverse() * ear_L_wp;
	const Eigen::Vector4f ear_R_in_mouth_coords = mouth_transform.inverse() * ear_R_wp;

	// calculate mouth/ear output angles and distances:
	getAzimuthAndElevationFromVector(mouth_in_ear_L_coords, spatialData.ear_L_to_mouth_azimuth, spatialData.ear_L_to_mouth_elevation);
	getAzimuthAndElevationFromVector(mouth_in_ear_R_coords, spatialData.ear_R_to_mouth_azimuth, spatialData.ear_R_to_mouth_elevation);

	getAzimuthAndElevationFromVector(ear_L_in_mouth_coords, spatialData.mouth_to_ear_L_azimuth, spatialData.mouth_to_ear_L_elevation);
	getAzimuthAndElevationFromVector(ear_R_in_mouth_coords, spatialData.mouth_to_ear_R_azimuth, spatialData.mouth_to_ear_R_elevation);

	spatialData.ear_L_to_mouth_dist = (mouth_wp - ear_L_wp).norm();
	spatialData.ear_R_to_mouth_dist = (mouth_wp - ear_R_wp).norm();

	std::cout << "ear_L_to_mouth_azimuth = \n" << spatialData.ear_L_to_mouth_azimuth << std::endl;
	std::cout << "ear_R_to_mouth_azimuth = \n" << spatialData.ear_R_to_mouth_azimuth << std::endl;
	std::cout << "mouth_to_ear_L_azimuth = \n" << spatialData.mouth_to_ear_L_azimuth << std::endl;
	std::cout << "mouth_to_ear_R_azimuth = \n" << spatialData.mouth_to_ear_R_azimuth << std::endl;


}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
