#pragma once

#include"LocalPara.h"
#include"SrcPara.h"
#include"SplineMap\Parafun.h"
#include"HybridMap.h"
#include"SharedData.h"
#include"SimplifyMesh.h"
class PolyMap
{
public:
	PolyMap();
	~PolyMap();

public:
	void load_mesh(const string& file_str);

	void init_with_partition();
	void init_with_tutte();

	void run_methods(int meth_type,int init_type);
	double calc_distortion(SharedData::DIS_WEIGHT dis_w);
	void update_pos_from_sD();

	void local_opt_with_fixed_boundary();
	bool local_opt_onebyone(vector<std::vector<int>>& phase,int max_iter_num_=20,double convergence_rate_=1e-4,bool is_fixed_b=false);
	void local_opt_1v(int vid);
	void local_decrease_1v(int vid, Eigen::Vector2d& p_g_);
	double local_energy_1v(int vid,const Eigen::Vector2d& p_);


	void add_scaffold_to_mesh();
	void delete_scaffold_from_mesh();


	double calc_cur_face_dis(vector<pair<double, int>>& e_dis_iso);
	void load_2d_mesh(const char* filename_);
	bool grow_from_seeds_upto_vnum(set<int>& fs_set, const set<int>& seeds_v, vector<int>& f_isaccessiable, int fnum_upper);
	void output_related_info(const string& file_);

	void strategy_our();

	void refinement_with_2duv(const char* filename_);

	void submesh_based_method();


private:

	template<typename T>
	void out_vec_file(const vector<T>& vec_, const string& file_);

public:
	SharedData sD;
	//Mesh mesh_whole;
	int V_N;
	int F_N;
	string model_name;
	//Eigen::VectorXd pos_of_mesh;

	Eigen::MatrixXd w_uv_backup;

	//solver related
	int VAR_N;
	int MAX_ITER_NUM = 5000;
	double convergence_rate = 1e-4;
	bool is_init_2d=false;

	//==================
	std::shared_ptr<Parafun> spline_solver = nullptr;
	std::shared_ptr<HybridMap> hybrid_solver = nullptr;
	std::shared_ptr<SrcPara> srcpara_solver = nullptr;

	bool is_spline_ = true;
};

template<typename T>
inline void PolyMap::out_vec_file(const vector<T>& vec_, const string & file_)
{
	ofstream ofs(file_, ios::trunc);
	for (const T &var : vec_)
		ofs << var << endl;
	ofs.close();
}
