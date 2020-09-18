#pragma once
#include"Common\Common.h"
#include"SplineMap\ScafData.h"
struct SharedData
{
public:
	SharedData();
public:
	enum DIS_WEIGHT
	{
		AREA, UNIFORM
	};

	DIS_WEIGHT dis_w_type = DIS_WEIGHT::AREA;


	//mesh
	string model_name;
	Mesh mesh;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
//	vector<vector<int>> VF;

	//vector<double> src_inv00;
	//vector<double> src_inv01;
	//vector<double> src_inv10;
	//vector<double> src_inv11;

	double total_area=1.0;
	double area_same_factor;

	//area-weight
	vector<double> src_area;
	vector<double> src_weight;

	//scaffold
	long mv_num, mf_num;
	long sv_num, sf_num;
	long V_num, F_num;
	//Eigen::MatrixXd m_V; // input initial mesh // mv_num*3
	//Eigen::MatrixXi m_T; // input initial mesh // mf_num*3

	Eigen::MatrixXd w_uv; // whole domain uv // (mv_num+sv_num)*2
	Eigen::MatrixXi s_T; // scaffold mesh faces // sf_num*3
	//Eigen::MatrixXi w_T; // whole mesh faces // (mf_num + sf_num)*3

	//Eigen::VectorXd m_M; // mesh area or volume
	//Eigen::VectorXd s_M; // scaffold area or volume
	//Eigen::VectorXd w_M; // area/volume weights for whole

	double mesh_measure; // whole area or volume
	int dim;
	double scaf_weight = 1.0;
	double area_threshold;

	double interval;

	Eigen::VectorXi frame_ids;
	Eigen::VectorXi internal_bnd;
	Eigen::MatrixXd rect_frame_V;

	// multi-chart support
	std::vector<int> component_sizes;//components num, faces num in every component
	std::vector<int> bnd_sizes;//loop num,vertex num in every loop 

	vector<Eigen::Vector4d> scaf_jac;


	//PP
	bool is_interp = false;
	vector<double> upd_inv00;
	vector<double> upd_inv01;
	vector<double> upd_inv10;
	vector<double> upd_inv11;
	double bound_distortion_K = 250.0;
	double intp_T_min;
	double changetocm_flag;

	////boundary
	std::set<int> order1_fs_b;
	std::vector<std::vector<int>> boundary_loops;//counterclockwise, boundary_outter
	vector<int> boundary_inner_loop;//boundary_inner
	vector<int> boundary_es_del;
	set<int> boundary_vs;
	vector<int> first_inner_bl;

	////opt
	double Energy_cur;


	////show
	vector<double> distortion;//show distortion & calc_order for retri
	vector<double> energy_record;

	vector<tuple<double, double, int>> energy_record3;

	//modify bad input
	double e_len_avg;
	double min_tri_area;

	//inner-seam for order1
//	std::set<int> order1_fs_innerseam;
	vector<std::vector<int>> phase;


	vector<set<int>> patches_fids;
	int patches_cur_id=0;

private:
	void init_V_F_VF();
	void rescale_and_init_area();
	void find_components(Eigen::VectorXi& F_flag);

public:
	void add_new_patch(const Eigen::MatrixXd&, const Eigen::MatrixXi&,
		const Eigen::RowVectorXd &center/*,const Eigen::VectorXd&*/);

	void mesh_improve(bool);
	void update_scaffold();
	void prepare_data(int loops_num);
	void init_src_inv();
	void restore_upd2src();
	void init_loops(int loops_num);
	void find_boundary_loops();
	void grow_from_seeds(set<int>& fs_set, const set<int>& seeds_v, int loops_num);
	void boundary_loop_for_part(set<int>& part, std::vector<int>& boundaryloop);
	void update_source_jacobian(const Eigen::VectorXd& pos_);
	void update_uv_from_outersrc(const Eigen::VectorXd& pos_);
	void adjust_scaf_weight(double new_weight);
	void handle_mintri();
	void update_scaf_jacobian();
	void prepare_data_with_partion(int loops_num);

	void write_obj_src_uv(const string &outstr);

	void local_opt_onebyone(vector<std::vector<int>>& phase, int max_iter_num_ = 20, double convergence_rate_ = 1e-4);
	void local_decrease_1v(int vid, Eigen::Vector2d & p_g_);
	void local_opt_1v(int vid);
	double local_energy_1v(int vid, const Eigen::Vector2d& p_);
	double calc_energy(const Eigen::MatrixXd & pos_cur);
	//droped
	void find_order1_region_near_boundary_face(set<int>& fs_set, int loops = 3);
	void find_order1_region_near_boundary_vertex(set<int>& fs_set, int loops = 3);
	void find_order1_region_near_boundary_by_radius(set<int>& fs_set, const set<int>& seeds_v, int loops_num=5);



};