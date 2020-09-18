#pragma once
#include"PardisoSolver.h"
#include<list>
#include <vector>
#include<set>
#include<fstream>
#include<math.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include"Eigen/Dense"
#include"Eigen/Sparse"
#include<time.h>
#include<omp.h>
#include<queue>
#include"MeshViewer\MeshDefinition.h"
#include<glog\logging.h>


using namespace std;
void orderTheBoundary(vector<list<int>>& order_boundary, const vector<int>& nextlocation);
bool growFromP(int p, set<int>& isused, list<int>& order_boundary, const vector<int>& nextlocation);

void map_vertices_to_circle(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd& UV);

void Tutte(int V_N, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd, const Eigen::MatrixXd& bnd_uv, Eigen::MatrixXd & uv_init);
void preCalc_pardiso(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, PardisoSolver & pardiso);

void writeObj(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const string& outfile);
void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryEdges);
void remove_unreferenced(int V_N, const Eigen::MatrixXi & Fsrc, Eigen::VectorXi& J, Eigen::MatrixXi& PF, Eigen::VectorXi& VI);

double get_smallest_pos_quad_zero(double a, double b, double c);
double newton_equation(const double & a, const double & b, const double & K);

void test_remove_unreferenced(const Eigen::MatrixXd & V, const Eigen::MatrixXi & Fsrc, Eigen::VectorXi& J, Eigen::MatrixXi& PF, Eigen::VectorXi& VI);

void graph_color(const Mesh& mesh, std::vector<std::vector<int>>& phase);

int preCalc_pardiso(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

