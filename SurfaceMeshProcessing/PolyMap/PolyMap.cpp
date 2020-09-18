#include "PolyMap.h"

PolyMap::PolyMap()
{
}

PolyMap::~PolyMap()
{

}

void PolyMap::load_mesh(const string & file_str)
{
	printf("loading mesh %s ......\n",file_str.c_str());
	OpenMesh::IO::read_mesh(sD.mesh, file_str);
	model_name = file_str;
	sD.model_name = file_str;

	cout << "load_mesh: " << file_str << " v_n: " << sD.mesh.n_vertices() << " f_n: " << sD.mesh.n_faces() << endl;
	LOG(INFO) << "load_mesh_name: " << file_str << " v_n: " << sD.mesh.n_vertices() << " f_n: " << sD.mesh.n_faces() << endl;
}

void PolyMap::init_with_partition()
{
	printf("init begin-----------(partition)\n");
	LOG(INFO) << "init begin-----------(partition)";
	clock_t t0 = clock();

	sD.find_boundary_loops();
	SimplifyMesh sply_mesh(sD);
	sply_mesh.color_src_mesh();
	sD.prepare_data_with_partion(5);

	V_N = sD.mv_num;
	F_N = sD.mf_num;
	sD.Energy_cur = calc_distortion(SharedData::DIS_WEIGHT::AREA);

	clock_t t1 = clock();
	double t_init = (double)(t1 - t0) / CLOCKS_PER_SEC;
	sD.energy_record.push_back(sD.Energy_cur);
	sD.energy_record3.emplace_back(sD.Energy_cur, t_init, 0);
	printf("init_with_partition finish-- %f --; energy_cur: %f ;\n", t_init, sD.Energy_cur);
	LOG(INFO) << "init_with_partition finish t-- " << t_init << " ; energy_cur: " << sD.Energy_cur;

}

void PolyMap::init_with_tutte()
{
	printf("init begin-----------(tutte)\n");
	LOG(INFO) << "init begin-----------(tutte)";
	clock_t t0 = clock();
	sD.find_boundary_loops();
	sD.prepare_data(5);

	V_N = sD.mv_num;
	F_N = sD.mf_num;
	sD.Energy_cur = calc_distortion(SharedData::DIS_WEIGHT::AREA);
	clock_t t1 = clock();
	double t_init = (double)(t1 - t0) / CLOCKS_PER_SEC;

	sD.energy_record.push_back(sD.Energy_cur);
	sD.energy_record3.emplace_back(sD.Energy_cur, t_init, 0);

	printf("init_with_tutte finish t: %f ; energy_cur: %f ;\n", t_init, sD.Energy_cur);
	LOG(INFO) << "init_with_tutte finish t: " << t_init << " ; energy_cur: " << sD.Energy_cur;

}

void PolyMap::run_methods(int meth_type, int init_type)
{
	if (init_type == 0)
	{
		init_with_tutte();
		writeObj(sD.w_uv, sD.F, "./" + model_name + "_tutte_para.obj");
		output_related_info("_tutte_src_uv");
	}
	else if (init_type == 1)
	{

		init_with_partition();
	}
	else
	{
		printf("Wrong init_type!\n");
		return;
	}

	if (meth_type == 0)
	{
		strategy_our();
	}
	else if (meth_type == 1|| meth_type == 2)
	{
		srcpara_solver.reset(new SrcPara(sD));
		srcpara_solver->set_method(meth_type - 1);
		srcpara_solver->run();
	}
	else
	{
		printf("Wrong meth_type!\n");
		return;
	}

	writeObj(sD.w_uv, sD.F, "./" + model_name + "_bi_para.obj");
	//output
	output_related_info("_bi_src_uv");

}

double PolyMap::calc_distortion(SharedData::DIS_WEIGHT dis_w)
{
//	auto& pos_of_mesh = hybrid_solver->pos_of_mesh;
	vector<vector<Eigen::RowVector2d>> pp;
	switch (dis_w)
	{
	case SharedData::AREA:
#pragma omp parallel for
		for (int i = 0; i < F_N; i++)
		{
			int f0 = sD.F(i, 0);
			int f1 = sD.F(i, 1);
			int f2 = sD.F(i, 2);

			//double q00 = pos_of_mesh[f0] - pos_of_mesh[f2];
			//double q01 = pos_of_mesh[f1] - pos_of_mesh[f2];
			//double q10 = pos_of_mesh[f0 + V_N] - pos_of_mesh[f2 + V_N];
			//double q11 = pos_of_mesh[f1 + V_N] - pos_of_mesh[f2 + V_N];

			double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
			double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
			double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
			double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);

			//double j00 = sD.src_inv00[i] * q00;
			//double j01 = sD.src_inv01[i] * q00 + sD.src_inv11[i] * q01;
			//double j10 = sD.src_inv00[i] * q10;
			//double j11 = sD.src_inv01[i] * q10 + sD.src_inv11[i] * q11;

			double j00 = sD.upd_inv00[i] * q00;
			double j01 = sD.upd_inv01[i] * q00 + sD.upd_inv11[i] * q01;
			double j10 = sD.upd_inv00[i] * q10;
			double j11 = sD.upd_inv01[i] * q10 + sD.upd_inv11[i] * q11;


			double det = j00 * j11 - j01 * j10;
			if (det <= 0)
			{
				cout << i << " lajilaji " << det << endl;
				LOG(INFO) << i << " lajilaji " << det << endl;
				pp.push_back(vector<Eigen::RowVector2d>{ sD.w_uv.row(f0), sD.w_uv.row(f1), sD.w_uv.row(f2) });

			}
			double E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			double E_2 = E_1 / (det*det);
			sD.distortion[i] = sD.src_area[i] * (E_1 + E_2);

			if (sD.distortion[i] > 1e4)
			{
				pp.push_back(vector<Eigen::RowVector2d>{ sD.w_uv.row(f0) , sD.w_uv.row(f1) ,sD.w_uv.row(f2) });
			}

		}
		break;
	case SharedData::UNIFORM:
#pragma omp parallel for
		for (int i = 0; i < F_N; i++)
		{
			int f0 = sD.F(i, 0);
			int f1 = sD.F(i, 1);
			int f2 = sD.F(i, 2);


			double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
			double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
			double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
			double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);


			//double q00 = pos_of_mesh[f0] - pos_of_mesh[f2];
			//double q01 = pos_of_mesh[f1] - pos_of_mesh[f2];
			//double q10 = pos_of_mesh[f0 + V_N] - pos_of_mesh[f2 + V_N];
			//double q11 = pos_of_mesh[f1 + V_N] - pos_of_mesh[f2 + V_N];

			double j00 = sD.upd_inv00[i] * q00;
			double j01 = sD.upd_inv01[i] * q00 + sD.upd_inv11[i] * q01;
			double j10 = sD.upd_inv00[i] * q10;
			double j11 = sD.upd_inv01[i] * q10 + sD.upd_inv11[i] * q11;

			double det = j00 * j11 - j01 * j10;

			double E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			double E_2 = E_1 / (det*det);
			sD.distortion[i] = E_1 + E_2;
		}
		break;
	default:
		break;
	}

//	out_vec_file<double>(sD.distortion, "dis.txt");


	ofstream ofs("./data/" + model_name + "x.txt", ios::trunc);

	for (size_t i = 0; i < pp.size(); i++)
	{
		for (size_t j = 0; j < pp[i].size(); j++)
		{
			ofs<<setprecision(10)<<fixed << pp[i][j] << " ";
		}
		ofs << endl;
	}

	ofs.close();

	double e_sum = 0.0;
	for (auto&var : sD.distortion)
		e_sum += var;
	return e_sum;
}

void PolyMap::update_pos_from_sD()
{
	sD.w_uv.conservativeResize(V_N, 2);
}

//0:accessiable/1:unaccessiable
bool PolyMap::grow_from_seeds_upto_vnum(set<int>& fs_set, const set<int>& seeds_v, vector<int>& f_isaccessiable, int fnum_upper)
{
	set<int> contain_vs;
	std::queue<int> Q;
	int actual_f_num = 0;
	for (auto&vid : seeds_v)
	{
		Q.emplace(vid);
		contain_vs.insert(vid);
		auto vid_h = sD.mesh.vertex_handle(vid);
		for (auto itvf = sD.mesh.vf_begin(vid_h); itvf != sD.mesh.vf_end(vid_h); itvf++)
		{
			int var = itvf->idx();
			if (var >= sD.mf_num)
			{
				fs_set.insert(var);
			}
			else if (f_isaccessiable[var] == 0)//==0
			{
				f_isaccessiable[var] = 1;
				fs_set.insert(var);
				actual_f_num++;
			}
			else if (f_isaccessiable[var] == 2)
			{
				fs_set.insert(var);
			}
		}

	}

	while (!Q.empty() && actual_f_num < fnum_upper)
	{
		auto var = Q.front();
		Q.pop();
		auto v_h = sD.mesh.vertex_handle(var);
		for (auto itvv = sD.mesh.vv_begin(v_h); itvv != sD.mesh.vv_end(v_h); itvv++)
		{
			int vid_ = itvv->idx();
			if (vid_ < sD.mv_num&&contain_vs.find(vid_) == contain_vs.end())
			{
				Q.emplace(vid_);
				contain_vs.insert(vid_);
				auto vid_h = sD.mesh.vertex_handle(vid_);
				for (auto itvf = sD.mesh.vf_begin(vid_h); itvf != sD.mesh.vf_end(vid_h); itvf++)
				{
					int vf_ = itvf->idx();
					if (vf_ >= sD.mf_num)
					{
						fs_set.insert(vf_);
					}
					else if (f_isaccessiable[vf_] == 0)//==0
					{
						f_isaccessiable[vf_] = 1;
						fs_set.insert(vf_);
						actual_f_num++;
					}
					else if (f_isaccessiable[vf_] == 2)
					{
						fs_set.insert(vf_);
					}
				}
			}
		}
	}



	if (fs_set.size() < 50)
		return false;

	double e_tmp_ = 0.0;
	double a_sum_ = 0.0;
	for (const auto &var_ : fs_set)
	{
		if (var_ < sD.mf_num)
		{
			e_tmp_ += sD.distortion[var_];
			a_sum_ += sD.src_weight[var_];
		}
	}
	e_tmp_ /= a_sum_;

	if (e_tmp_ < sD.Energy_cur - 0.2)
		return false;
	else
		return true;
	

}

void PolyMap::output_related_info(const string & file_)
{
	ofstream ofs("./" + model_name + "_plot.txt", ios::trunc);
	for (auto &var : sD.energy_record3)
		ofs << get<0>(var) << " " << get<1>(var) << " " << get<2>(var) << endl;
	ofs.close();

	sD.w_uv.conservativeResize(sD.mv_num, 2);
	{
		Eigen::MatrixXd& V_in = sD.w_uv;
		Eigen::RowVector2d f_bb_max = V_in.colwise().maxCoeff();
		Eigen::RowVector2d f_bb_min = V_in.colwise().minCoeff();

		auto bbox = f_bb_max - f_bb_min;
		double max_size = bbox.maxCoeff();
		for (int i = 0; i < V_in.rows(); i++)
		{
			Eigen::RowVector2d v_p = V_in.row(i);
			v_p = v_p - f_bb_min;
			v_p /= max_size;
			V_in.row(i) = v_p;
		}
	}
	printf("----------------write para_with_UV-----------\n");
	sD.write_obj_src_uv("./" + model_name + file_+".obj");

}

void PolyMap::strategy_our()
{
	printf("strategy_our begin-----------\n");
	LOG(INFO) << "strategy_our begin-----------";
	clock_t t0 = clock();
	spline_solver = make_shared<Parafun>(sD);
	hybrid_solver = make_shared<HybridMap>(sD);
	hybrid_solver->prepare_data();
	if (sD.Energy_cur > 1e5)
	{
		vector<std::vector<int>> phase;
		graph_color(sD.mesh, phase);
		bool is_conv_1by1_ = false;
		while (sD.Energy_cur > 1e4&&!is_conv_1by1_)
		{
			is_conv_1by1_=local_opt_onebyone(phase,10,1e-4,true);
		}
	}
	else
	{
		vector<std::vector<int>> phase;
		graph_color(sD.mesh, phase);
		local_opt_onebyone(phase, 10, 1e-4, true);
	}
	get<2>(sD.energy_record3.back()) = 1;

	clock_t t1 = clock();
	double t_pw_ = (double)(t1 - t0) / CLOCKS_PER_SEC;


	spline_solver->iteratively_BPE_spline();
	update_pos_from_sD();
	get<2>(sD.energy_record3.back()) = 1;
	clock_t t2 = clock();
	double t_sp_ = (double)(t2 - t1) / CLOCKS_PER_SEC;


	hybrid_solver->autorun(1e-4);
	get<2>(sD.energy_record3.back()) = 1;
	clock_t t3 = clock();
	double t_hy_ = (double)(t3 - t2) / CLOCKS_PER_SEC;


	vector<std::vector<int>> phase;
	graph_color(sD.mesh, phase);
	local_opt_onebyone(phase, 20);
	phase.clear();
	local_opt_with_fixed_boundary();
	clock_t t4 = clock();
	double t_loc_ = (t4 - t3) / CLOCKS_PER_SEC;

	cout << "t_pw: " << t_pw_ << " t_sp: " << t_sp_ << " t_hy: " << t_hy_ << " t_loc: " << t_loc_ << endl;
	LOG(INFO) << "t_pw: " << t_pw_ << " t_sp: " << t_sp_ << " t_hy: " << t_hy_ << " t_loc: " << t_loc_;

}

void PolyMap::refinement_with_2duv(const char * filename_)
{
	printf("refinement_with_2duv--begin-----\n");

	sD.find_boundary_loops();
	sD.prepare_data_with_partion(5);

	V_N = sD.mv_num;
	F_N = sD.mf_num;
	load_2d_mesh(filename_);

	sD.Energy_cur = calc_distortion(SharedData::DIS_WEIGHT::AREA);

	printf("prepare_data finish-- energy_cur: %f ;\n", sD.Energy_cur);


	vector<std::vector<int>> phase;
	graph_color(sD.mesh, phase);
	//local_opt_onebyone(phase,100000,1e-6);
	local_opt_onebyone(phase, 20);
	local_opt_with_fixed_boundary();

	sD.w_uv.conservativeResize(sD.mv_num, 2);
	writeObj(sD.w_uv, sD.F, "./data/" + model_name+"_localopt.obj");
}

void PolyMap::submesh_based_method()
{
	clock_t t1, t2, t3, tc0, tc1;
	printf("prepare_data begin-----------(origin_tutte)\n");
	LOG(INFO) << "submesh_based_method prepare_data begin-----------(origin_tutte)";
	t1 = clock();
	sD.find_boundary_loops();
	sD.prepare_data(5);

	V_N = sD.mv_num;
	F_N = sD.mf_num;
	t2 = clock();

	sD.Energy_cur = calc_distortion(SharedData::DIS_WEIGHT::AREA);
	sD.energy_record.push_back(sD.Energy_cur);
	sD.energy_record3.emplace_back(sD.Energy_cur, 0, 0);
	printf("prepare_data finish-- %f --; energy_cur: %f ;\n", (t2 - t1) / 1000.0, sD.Energy_cur);
	LOG(INFO) << "prepare_data finish t-- " << (t2 - t1) / 1000.0 << " ; energy_cur: " << sD.Energy_cur;

	local_opt_with_fixed_boundary();

	{
		ofstream ofs("./data/" + model_name + "_plot.txt", ios::trunc);
		for (auto &var : sD.energy_record3)
			ofs << get<0>(var) << " " << get<1>(var) << " " << get<2>(var) << endl;
		ofs.close();

	}

	out_vec_file<double>(sD.energy_record, "./data/" + model_name + "_energy_submesh.txt");
	sD.w_uv.conservativeResize(sD.mv_num, 2);
	{
		Eigen::MatrixXd& V_in = sD.w_uv;
		Eigen::RowVector2d f_bb_max = V_in.colwise().maxCoeff();
		Eigen::RowVector2d f_bb_min = V_in.colwise().minCoeff();

		auto bbox = f_bb_max - f_bb_min;
		double max_size = bbox.maxCoeff();
		LOG(INFO) << "if want restore, then times the factor: " << max_size;
		std::cout << max_size << std::endl;
		std::cout << f_bb_max << std::endl;
		std::cout << f_bb_min << std::endl;
		for (int i = 0; i < V_in.rows(); i++)
		{
			Eigen::RowVector2d v_p = V_in.row(i);
			v_p = v_p - f_bb_min;
			v_p /= max_size;
			V_in.row(i) = v_p;
		}
	}
	printf("----------------write runsrc_with_UV-----------\n");
	sD.write_obj_src_uv("./data/" + model_name + "_submesh_with_UV.obj");

}

void PolyMap::local_opt_with_fixed_boundary()
{
	vector<pair<double, int>> e_dis_iso(F_N);
	sD.Energy_cur = calc_cur_face_dis(e_dis_iso);
	printf("Enter the local refinement phase(local_opt_with_fixed_boundary) with init energy: %f.............\n", sD.Energy_cur);

	LOG(INFO) << "Enter the local refinement phase(local_opt_with_fixed_boundary) with init energy: " << sD.Energy_cur;

	int local_size_ = max(1000000, min(F_N / 20, 1000000));

	double e_pre_;
	double cur_conv_rate = 1.0;
	for (size_t loc_refine_ = 0; loc_refine_ < 3; loc_refine_++)
	{
		clock_t t1 = clock();

		sD.mesh_improve(true);
		sD.adjust_scaf_weight(sD.Energy_cur / sD.sf_num / 100.0);
		sD.update_scaf_jacobian();
		add_scaffold_to_mesh();

		e_pre_ = sD.Energy_cur;
		vector<int> f_visit(F_N, 0);

//#pragma omp parallel for
		for (int i = 0; i < F_N; i++)
		{
			if (f_visit[e_dis_iso[i].second] != 0)//==1
				continue;
			set<int> seeds_v;
			seeds_v.insert(sD.F(e_dis_iso[i].second, 0));
			seeds_v.insert(sD.F(e_dis_iso[i].second, 1));
			seeds_v.insert(sD.F(e_dis_iso[i].second, 2));
			set<int> Fs_p;
			bool is_executed_;
//#pragma omp critical
			{
				is_executed_ = grow_from_seeds_upto_vnum(Fs_p, seeds_v, f_visit, local_size_);
				//int f_num_unvisited = 0;
				//for (const auto& var : f_visit)
				//{
				//	if (var > 0)
				//		f_num_unvisited++;
				//}
				//cout << "f_num_unvisited " << f_num_unvisited << endl;
			}
			if (is_executed_)
			{
				LocalPara lp_(sD, Fs_p, f_visit);
				lp_.run();
			}
		}
		delete_scaffold_from_mesh();
		clock_t t2 = clock();
		sD.Energy_cur = calc_cur_face_dis(e_dis_iso);
		sD.energy_record.push_back(sD.Energy_cur);
		sD.energy_record3.emplace_back(sD.Energy_cur, (t2 - t1) / 1000.0, 0);
		cur_conv_rate = abs(e_pre_ - sD.Energy_cur) / e_pre_;
		
		printf("Local-Refinement iterations: %d ,Time: %f ,Energy: %f \n", loc_refine_,(t2-t1)/1000.0, sD.Energy_cur);
		LOG(INFO) << "Local-Refinement iterations: " << loc_refine_ << " ,Time: " << (t2 - t1) / 1000.0 << " ,Energy: " << sD.Energy_cur;
		if (cur_conv_rate < convergence_rate)
			break;
		if (loc_refine_>0&&loc_refine_ % 50 == 0)
		{
			sD.w_uv.conservativeResize(sD.mv_num, 2);
			{
				Eigen::MatrixXd& V_in = sD.w_uv;
				Eigen::RowVector2d f_bb_max = V_in.colwise().maxCoeff();
				Eigen::RowVector2d f_bb_min = V_in.colwise().minCoeff();
				auto bbox = f_bb_max - f_bb_min;
				double max_size = bbox.maxCoeff();
				LOG(INFO) << "if want restore, then times the factor: " << max_size;
				std::cout << max_size << std::endl;
				std::cout << f_bb_max << std::endl;
				std::cout << f_bb_min << std::endl;
				for (int i = 0; i < V_in.rows(); i++)
				{
					Eigen::RowVector2d v_p = V_in.row(i);
					v_p = v_p - f_bb_min;
					v_p /= max_size;
					V_in.row(i) = v_p;
				}
			}
			printf("----------------write runsrc_with_UV-----------\n");
			sD.write_obj_src_uv("./data/" + model_name + to_string(loc_refine_)+"_submesh_with_UV.obj");
		}

	}

}

bool PolyMap::local_opt_onebyone(vector<std::vector<int>>& phase, int max_iter_num_, double convergence_rate_, bool is_fixed_b)
{
	vector<pair<double, int>> e_dis_iso(F_N);
	sD.Energy_cur = calc_cur_face_dis(e_dis_iso);

	double e_pre_;
	double cur_conv_rate = 1.0;

	if (is_fixed_b)
	{
		printf("Enter the local refinement phase(onebyone_fixed_boundary) with init energy: %f.............\n", sD.Energy_cur);
		LOG(INFO) << "Enter the local refinement phase(onebyone_fixed_boundary) with init energy: " << sD.Energy_cur;
		for (size_t loc_refine_ = 0; loc_refine_ < max_iter_num_; loc_refine_++)
		{
			clock_t t1 = clock();
			e_pre_ = sD.Energy_cur;
			for (const auto&c_vids : phase)
			{
#pragma omp parallel for
				for (int i = 0; i < c_vids.size(); i++)
				{
					if(sD.boundary_vs.find(c_vids[i])==sD.boundary_vs.end())
						local_opt_1v(c_vids[i]);
				}
			}
			clock_t t2 = clock();
			sD.Energy_cur = calc_cur_face_dis(e_dis_iso);
			sD.energy_record.push_back(sD.Energy_cur);
			sD.energy_record3.emplace_back(sD.Energy_cur, (t2 - t1) / 1000.0, 0);
			cur_conv_rate = (e_pre_ - sD.Energy_cur) / e_pre_;
			
			printf("Local-Refinement iterations: %d ,Time: %f ,Energy: %f \n", loc_refine_, (t2 - t1) / 1000.0, sD.Energy_cur);
			LOG(INFO) << "Local-Refinement iterations: " << loc_refine_ << " Time: " << (t2 - t1) / 1000.0 << " Energy: " << sD.Energy_cur;
			if (cur_conv_rate < convergence_rate_)
			{
				return true;
			}
		}
	}
	else
	{
		printf("Enter the local refinement phase(onebyone_free_boundary) with init energy: %f.............\n", sD.Energy_cur);
		LOG(INFO) << "Enter the local refinement phase(onebyone_free_boundary) with init energy: " << sD.Energy_cur;
		for (size_t loc_refine_ = 0; loc_refine_ < max_iter_num_; loc_refine_++)
		{
			clock_t t1 = clock();

			sD.mesh_improve(true);
			sD.adjust_scaf_weight(sD.Energy_cur / sD.sf_num / 100.0);
			sD.update_scaf_jacobian();
			add_scaffold_to_mesh();

			e_pre_ = sD.Energy_cur;

			for (const auto&c_vids : phase)
			{
#pragma omp parallel for
				for (int i = 0; i < c_vids.size(); i++)
				{
					local_opt_1v(c_vids[i]);
				}
			}
			delete_scaffold_from_mesh();
			clock_t t2 = clock();
			sD.Energy_cur = calc_cur_face_dis(e_dis_iso);
			sD.energy_record.push_back(sD.Energy_cur);
			sD.energy_record3.emplace_back(sD.Energy_cur, (t2 - t1) / 1000.0, 0);
			cur_conv_rate = (e_pre_ - sD.Energy_cur) / e_pre_;
			printf("Local-Refinement iterations: %d ,Time: %f ,Energy: %f \n", loc_refine_, (t2 - t1) / 1000.0, sD.Energy_cur);
			LOG(INFO) << "Local-Refinement iterations: " << loc_refine_ << " Time: " << (t2 - t1) / 1000.0 << " Energy: " << sD.Energy_cur;
			if (cur_conv_rate < convergence_rate_)
			{
				return true;
			}
		}
	}
	return false;
}

double PolyMap::local_energy_1v(int vid, const Eigen::Vector2d& p_)
{
	double e_sum = 0.;
	auto seed_h = sD.mesh.vertex_handle(vid);
	for (auto itvf = sD.mesh.vf_begin(seed_h); itvf != sD.mesh.vf_end(seed_h); itvf++)
	{
		int fid_ = itvf->idx();
		int f0, f1, f2;
		double p00, p01, p10, p11, area_now;
		if (fid_ < sD.mf_num)
		{
			f0 = sD.F(fid_, 0);
			f1 = sD.F(fid_, 1);
			f2 = sD.F(fid_, 2);
			p00 = sD.upd_inv00[fid_];
			p01 = sD.upd_inv01[fid_];
			p10 = sD.upd_inv10[fid_];
			p11 = sD.upd_inv11[fid_];
			area_now = sD.src_weight[fid_];
		}
		else
		{
			//continue;
			int fid_i = fid_ - sD.mf_num;
			f0 = sD.s_T(fid_i, 0);
			f1 = sD.s_T(fid_i, 1);
			f2 = sD.s_T(fid_i, 2);
			p00 = sD.scaf_jac[fid_i][0];
			p01 = sD.scaf_jac[fid_i][1];
			p10 = sD.scaf_jac[fid_i][2];
			p11 = sD.scaf_jac[fid_i][3];
			area_now = sD.scaf_weight;
		}

		double q00, q01, q10, q11;
		if (vid == f0)
		{
			q00 = p_(0) - sD.w_uv(f2, 0);
			q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
			q10 = p_(1) - sD.w_uv(f2, 1);
			q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);
		}
		else if (vid == f1)
		{
			q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
			q01 = p_(0) - sD.w_uv(f2, 0);
			q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
			q11 = p_(1) - sD.w_uv(f2, 1);
		}
		else
		{
			q00 = sD.w_uv(f0, 0) - p_(0);
			q01 = sD.w_uv(f1, 0) - p_(0);
			q10 = sD.w_uv(f0, 1) - p_(1);
			q11 = sD.w_uv(f1, 1) - p_(1);
		}

		double j00 = p00 * q00 + p10 * q01;
		double j01 = p01 * q00 + p11 * q01;
		double j10 = p00 * q10 + p10 * q11;
		double j11 = p01 * q10 + p11 * q11;
		double det = j00 * j11 - j01 * j10;
		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
		double e_ = area_now * (tr + tr / (det*det));

		e_sum += e_;
	}
	return e_sum;
}

void PolyMap::local_opt_1v(int vid)
{
	Eigen::Vector2d p_g_;
	local_decrease_1v(vid, p_g_);
	auto seed_h = sD.mesh.vertex_handle(vid);
	//max-step
	double max_step_ = numeric_limits<double>::infinity();
	for (auto itvoh = sD.mesh.voh_begin(seed_h); itvoh != sD.mesh.voh_end(seed_h); itvoh++)
	{
		int f1 = sD.mesh.to_vertex_handle(*itvoh).idx();
		int f2 = sD.mesh.to_vertex_handle(sD.mesh.next_halfedge_handle(*itvoh)).idx();
		double q00 = sD.w_uv(vid, 0) - sD.w_uv(f2, 0);
		double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		double q10 = sD.w_uv(vid, 1) - sD.w_uv(f2, 1);
		double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);

		double b = q00 * q11 - q10 * q01;
		double a = p_g_(0)*q11 - p_g_(1)*q01;
		double t_tmp = -b / a;
		if (t_tmp>0&&max_step_ > t_tmp)
		{
			max_step_ = t_tmp;
		}
	}
	double step_ = 0.9*max_step_;
	//line-search

	double h_ = 0.5;
	double c_ = 0.2;
	double tt = -(p_g_.squaredNorm());
	Eigen::Vector2d p_1v_(sD.w_uv(vid, 0), sD.w_uv(vid, 1));
	double ex = local_energy_1v(vid, p_1v_);
//	printf("ex==: %f\n", ex);

	Eigen::Vector2d pos_new = p_1v_ + step_ * p_g_;
	double e = local_energy_1v(vid, pos_new);
//	printf("e==: %f\n", e);

	double e_cache = e;
	double step_cache = step_;
	int count_c = 0;
	while (e > ex + step_ * c_*tt&&count_c<5)
	{
		step_ *= h_;
		pos_new = p_1v_ + step_ * p_g_;
		e = local_energy_1v(vid, pos_new);
		count_c++;
		if (e < e_cache)
		{
			e_cache = e;
			step_cache = step_;
		}
	//	printf("e==: %f\n", e);
	}
	if (e_cache < e)
	{
		step_ = step_cache;
	}

	sD.w_uv(vid, 0) = sD.w_uv(vid, 0) + step_ * p_g_(0);
	sD.w_uv(vid, 1) = sD.w_uv(vid, 1) + step_ * p_g_(1);

}

//gradient decrease
void PolyMap::local_decrease_1v(int vid, Eigen::Vector2d & p_g_)
{
	p_g_.setZero();
	auto seed_h = sD.mesh.vertex_handle(vid);
	for (auto itvf = sD.mesh.vf_begin(seed_h); itvf != sD.mesh.vf_end(seed_h); itvf++)
	{
		int fid_ = itvf->idx();
		int f0, f1, f2;
		double p00, p01, p10, p11, area_now;
		if (fid_ < sD.mf_num)
		{
			f0 = sD.F(fid_, 0);
			f1 = sD.F(fid_, 1);
			f2 = sD.F(fid_, 2);
			p00 = sD.upd_inv00[fid_];
			p01 = sD.upd_inv01[fid_];
			p10 = sD.upd_inv10[fid_];
			p11 = sD.upd_inv11[fid_];
			area_now = sD.src_weight[fid_];
		}
		else
		{
	//		continue;
			int fid_i = fid_ - sD.mf_num;
			f0 = sD.s_T(fid_i, 0);
			f1 = sD.s_T(fid_i, 1);
			f2 = sD.s_T(fid_i, 2);
			p00 = sD.scaf_jac[fid_i][0];
			p01 = sD.scaf_jac[fid_i][1];
			p10 = sD.scaf_jac[fid_i][2];
			p11 = sD.scaf_jac[fid_i][3];
			area_now = sD.scaf_weight;
		}
		double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
		double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
		double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);

		double j00 = p00 * q00 + p10 * q01;
		double j01 = p01 * q00 + p11 * q01;
		double j10 = p00 * q10 + p10 * q11;
		double j11 = p01 * q10 + p11 * q11;
		double det = j00 * j11 - j01 * j10;
		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;

		double det_2plus1 = 2.0 / (det * det) + 2.0;
		double det_3tr2 = 2.0*tr / (det * det * det);
		double r0 = area_now * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = area_now * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = area_now * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = area_now * (j11 * det_2plus1 - det_3tr2 * j00);

		if (vid == f0)
		{
			p_g_(0) -= (r0*p00 + r1 * p01);
			p_g_(1) -= (r2*p00 + r3 * p01);
		}
		else if (vid == f1)
		{
			p_g_(0) -= (r0*p10 + r1 * p11);
			p_g_(1) -= (r2*p10 + r3 * p11);
		}
		else
		{
			p_g_(0) += (r0*(p00 + p10) + r1 * (p01 + p11));
			p_g_(1) += (r2*(p00 + p10) + r3 * (p01 + p11));
		}
	}
}

void PolyMap::add_scaffold_to_mesh()
{
	for (size_t i = sD.mv_num; i < sD.V_num; i++)
	{
		OpenMesh::Vec3d pos_(sD.w_uv(i, 0), sD.w_uv(i, 1), 0.0);
		sD.mesh.add_vertex(pos_);
	}
	for (size_t i = 0; i < sD.sf_num; i++)
	{
		std::vector<OpenMesh::VertexHandle> face_v_hs;
		for (size_t j = 0; j < sD.s_T.cols(); j++)
		{
			face_v_hs.push_back(sD.mesh.vertex_handle(sD.s_T(i, j)));
		}
		sD.mesh.add_face(face_v_hs);
	}
	 
	cout << "Local= after add scaffold to mesh, V_N: " << sD.mesh.n_vertices() << " F_N: " << sD.mesh.n_faces() << " ( "<<sD.V_num<<" , "<<sD.F_num<<" );" << endl;

}

void PolyMap::delete_scaffold_from_mesh()
{
	for (size_t i = sD.mf_num; i < sD.F_num; i++)
	{
		//sD.mesh.delete_vertex(sD.mesh.vertex_handle(i));
		sD.mesh.delete_face(sD.mesh.face_handle(i));
	}
	sD.mesh.garbage_collection();
	cout << "Local= after delete scaffold to mesh, V_N: " << sD.mesh.n_vertices() << " F_N: " << sD.mesh.n_faces() << " ( " << sD.mv_num << " , " << sD.mf_num << " );" << endl;
}

double PolyMap::calc_cur_face_dis(vector<pair<double, int>>& e_dis_iso)
{
//	auto& pos_of_mesh = hybrid_solver->pos_of_mesh;
	double e_sum_ = 0.0;
#pragma omp parallel for reduction(+:e_sum_)
	for (int i = 0; i < F_N; i++)
	{
		int f0 = sD.F(i, 0);
		int f1 = sD.F(i, 1);
		int f2 = sD.F(i, 2);
		double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
		double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
		double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);

		double p00 = sD.upd_inv00[i];
		double p01 = sD.upd_inv01[i];
		double p10 = sD.upd_inv10[i];
		double p11 = sD.upd_inv11[i];
		double j00 = p00 * q00 + p10 * q01;
		double j01 = p01 * q00 + p11 * q01;
		double j10 = p00 * q10 + p10 * q11;
		double j11 = p01 * q10 + p11 * q11;
		double det = j00 * j11 - j01 * j10;
		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
		e_dis_iso[i].first = sD.src_weight[i]*(tr + tr / (det*det));
		e_dis_iso[i].second = i;
		sD.distortion[i]= e_dis_iso[i].first;
		e_sum_ += e_dis_iso[i].first;
	}
	sort(e_dis_iso.begin(), e_dis_iso.end(), greater<pair<double, int>>());
	return e_sum_;
}

void PolyMap::load_2d_mesh(const char * filename_)
{
	Mesh mesh_2d_;
	OpenMesh::IO::read_mesh(mesh_2d_, filename_);
	int v_2d_N = mesh_2d_.n_vertices();
	int f_2d_N = mesh_2d_.n_faces();
	sD.w_uv.resize(v_2d_N, 2);
	for (auto itv = mesh_2d_.vertices_begin(); itv != mesh_2d_.vertices_end(); itv++)
	{
		auto pos_ = mesh_2d_.point(*itv);
		sD.w_uv(itv->idx(), 0) = pos_[0];
		sD.w_uv(itv->idx(), 1) = pos_[1];
	}
	is_init_2d = true;
	printf("Load-2d-Mesh.......................V: %d;F: %d;\n", v_2d_N, f_2d_N);
}
