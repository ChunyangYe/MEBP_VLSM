#pragma once
#include<Eigen/Core>
#include<triangle.h>
void triangulate(const Eigen::MatrixXd &bnd_pts,const double& area_threshold,Eigen::MatrixXd & pts,Eigen::MatrixXi & FV);
void triangulate(const Eigen::MatrixXd &bnd_pts, const Eigen::MatrixXi &E, const Eigen::MatrixXd &hole, const double & area_threshold, Eigen::MatrixXd & pts, Eigen::MatrixXi & FV);
void triangulate(const Eigen::MatrixXd &bnd_pts, const Eigen::MatrixXi &E, const Eigen::MatrixXd &hole, Eigen::MatrixXd & pts, Eigen::MatrixXi & FV);


