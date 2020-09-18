#include "HybridMap.h"

HybridMap::HybridMap(SharedData& sharedd_):sD(sharedd_)
{
	
}

HybridMap::~HybridMap()
{

}

void HybridMap::solve_equation()
{
	clock_t t1, t2,t_beg,t_end;
	t_beg = clock();
	t1 = clock();
	prepare_pardiso();
	//prepare_pardiso_directly();//nonono
	t2 = clock();
	double t_pre = (t2 - t1) / 1000.0;
	//
	if (is_slim_meth)
	{
		assembly_SLIM_para();
		//assembly_directly_SLIM_para();//nonono
	}
	else
	{
		//assembly_hessian_grad_plustype();
		assembly_CM_para();
		//assembly_directly_CM_para();//nonono
	}

	t1 = clock();
	double t_ass = (t1 - t2) / 1000.0;

	pardiso_solver->pardiso_init();
	t2 = clock();
	double t_init = (t2 - t1) / 1000.0;

	pardiso_solver->factorize();
	t1 = clock();
	double t_factorize = (t1 - t2) / 1000.0;
	pardiso_solver->pardiso_solver();
	t2 = clock();
	double t_solve = (t2 - t1) / 1000.0;

	if (is_slim_meth)
	{
		printf("SLIM==");
		LOG(INFO) << "SLIM== ";
	}
	else
	{
		printf("CM==");
		LOG(INFO) << "CM== ";
	}
	printf("Pardiso prepare:%f //samplePoints: %d assembly:%f //solve equation;%f = init:%f + factorize:%f + pardiso_solve:%f ;\n", t_pre, (F_N - fix_F_N)*sample_num_per_tri +order1_fix_F_N-fix_F_N, t_ass, (t_init+ t_factorize+t_solve),t_init, t_factorize, t_solve);
	LOG(INFO) << "Pardiso prepare: "<< t_pre <<" //samplePoints: "<< (F_N - fix_F_N)*sample_num_per_tri + order1_fix_F_N - fix_F_N <<" assembly: "<< t_ass <<" //solve equation; "<< (t_init + t_factorize + t_solve) <<" = init: "<< t_init <<" + factorize: "<< t_factorize <<" + pardiso_solve: "<< t_solve <<" ;";
	vector<double>& ctrl_d_direction = pardiso_solver->result;
	Eigen::MatrixXd p_d_,p_g_;//pAs*2

	vector<double> t_vecc;
	t_vecc.emplace_back(clock() / 1000.0);
	calc_src_decrease_direction(ctrl_d_direction, p_d_);
	calc_src_gradient(pardiso_solver->rhs, p_g_);
	t_vecc.emplace_back(clock() / 1000.0);

	double max_step_ = max_step(p_d_);
	t_vecc.emplace_back(clock() / 1000.0);

	double step_;
	if(is_slim_meth)
		step_ = min(1.0, 0.8 * max_step_);
	else
		step_ = 0.95*max_step_;

	t_vecc.emplace_back(clock() / 1000.0);

	backtracking_line_search(p_d_, p_g_, step_);
	t_vecc.emplace_back(clock() / 1000.0);

	calc_energy(sD.w_uv);

	sD.w_uv += step_ * p_d_.topRows(V_N);
	for (size_t i = 0; i < 2*VAR_N; i++)
	{
		variables_[i] += step_ * ctrl_d_direction[i];
	}	
	sD.Energy_cur = calc_energy(sD.w_uv);
	t_vecc.emplace_back(clock() / 1000.0);

	t_end = clock();

	printf("calc_src_decrease_direction//max_step//calc_src_gradient//backtracking_line_search//calc-energy \n");
	for (int ti = 0; ti < t_vecc.size() - 1; ti++)
	{
		cout << t_vecc[ti + 1] - t_vecc[ti] << "//";
	}
	cout << endl;

	printf("Iteration finish(iter_time %f ;Max step: %f , Step length: %f , Cur_Energy: %f );\n",(t_end-t_beg)/1000.0, max_step_, step_, sD.Energy_cur);
	LOG(INFO)<<"Iteration finish(iter_time " << (t_end - t_beg) / 1000.0 << " ;Max step: "<< max_step_ <<" , Step length: "<< step_ <<"  , Cur_Energy: "<< sD.Energy_cur <<"  );";

}

void HybridMap::calc_src_decrease_direction(const vector<double>& c_d_, Eigen::MatrixXd& p_d_)
{
	pAs = V_N + order1_fix_V_N - fix_V_N;
	p_d_.setZero(pAs, 2);
#pragma omp parallel for
	for (int i = 0; i < V_N; i++)
	{
		int s_vid = fix_v_src2s[i];
		if (s_vid != -1)
		{
			p_d_(i,0) = c_d_[s_vid];
			p_d_(i,1) = c_d_[s_vid + VAR_N];
		}
		else
		{
			auto& fid_pos = fine2coarse_v_id_pos[i];
			auto& ids_bases = s_mesh.data(s_mesh.face_handle(fid_pos.first)).var_ids_bases;
			double du = 0.0;
			double dv = 0.0;
			for (size_t j = 0; j < ids_bases.size(); j++)
			{
				auto& res_ = ids_bases[j].second(fid_pos.second[0], fid_pos.second[1]);
				du += std::get<0>(res_)*c_d_[ids_bases[j].first];
				dv += std::get<0>(res_)*c_d_[ids_bases[j].first + VAR_N];
			}
			p_d_(i, 0) = du;
			p_d_(i, 1) = dv;
		}
	}

	for (int i = V_N; i < pAs; i++)
	{
		int sm_vid_ = i - V_N + fix_V_N;
		if (mid2varid[sm_vid_] != -1)
		{
			p_d_(i, 0) = c_d_[mid2varid[sm_vid_]];
			p_d_(i, 1) = c_d_[mid2varid[sm_vid_] + VAR_N];
		}
	}

	scaf_F.resize(order1_fix_F_N - fix_F_N,3);
	for (int i = fix_F_N; i < order1_fix_F_N; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (s_F(i, j) < fix_V_N)
				scaf_F(i - fix_F_N, j) = fix_v_s2src[s_F(i, j)];
			else
				scaf_F(i - fix_F_N, j) = s_F(i, j) - fix_V_N + V_N;
		}
	}

	////test
	//vector<double> ppd(2 * V_N);
	//for (size_t i = 0; i < 2*V_N; i++)
	//{
	//	ppd[i] = p_d_[i];
	//}
	//out_vec_file<double>(ppd, "d_dir.txt");

}

void HybridMap::calc_src_gradient(const vector<double>& c_g_, Eigen::MatrixXd & p_g_)
{
#if 0
	p_g_.resize(V_N * 2);
	for (int i = 0; i < V_N; i++)
	{
		int s_vid = fix_v_src2s[i];
		if (s_vid != -1)
		{
			p_g_[i] = c_g_[s_vid];
			p_g_[i + V_N] = c_g_[s_vid + VAR_N];
		}
		else
		{
			auto& fid_pos = fine2coarse_v_id_pos[i];
			auto& ids_bases = s_mesh.data(s_mesh.face_handle(fid_pos.first)).var_ids_bases;
			double du = 0.0;
			double dv = 0.0;
			for (size_t j = 0; j < ids_bases.size(); j++)
			{
				auto& res_ = ids_bases[j].second(fid_pos.second[0], fid_pos.second[1]);
				du += std::get<0>(res_)*c_g_[ids_bases[j].first];
				dv += std::get<0>(res_)*c_g_[ids_bases[j].first + VAR_N];
			}
			p_g_[i] = du;
			p_g_[i + V_N] = dv;
		}
	}

#else
	p_g_.setZero(pAs, 2);

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

		double det_2plus1 = 2.0 / det / det + 2.0;
		double det_3tr2 = 2.0*tr / det / det / det;
		double r0 = sD.src_weight[i] * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = sD.src_weight[i] * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = sD.src_weight[i] * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = sD.src_weight[i] * (j11 * det_2plus1 - det_3tr2 * j00);

		p_g_(f0,0) -= (r0*p00 + r1 * p01);
		p_g_(f1,0) -= (r0*p10 + r1 * p11);
		p_g_(f2,0) += (r0*(p00 + p10) + r1 * (p01 + p11));
		p_g_(f0,1) -= (r2*p00 + r3 * p01);
		p_g_(f1,1) -= (r2*p10 + r3 * p11);
		p_g_(f2,1) += (r2*(p00 + p10) + r3 * (p01 + p11));
	}

	for (int i = fix_F_N; i < order1_fix_F_N; i++)
	{
		int i_p = i - fix_F_N;
		int f0 = s_F(i, 0);
		int f1 = s_F(i, 1);
		int f2 = s_F(i, 2);
		double q00 = s_V(f0, 0) - s_V(f2, 0);
		double q01 = s_V(f1, 0) - s_V(f2, 0);
		double q10 = s_V(f0, 1) - s_V(f2, 1);
		double q11 = s_V(f1, 1) - s_V(f2, 1);
		double p00 = sD.scaf_jac[i_p][0];
		double p01 = sD.scaf_jac[i_p][1];
		double p10 = sD.scaf_jac[i_p][2];
		double p11 = sD.scaf_jac[i_p][3];
		double j00 = p00 * q00 + p10 * q01;
		double j01 = p01 * q00 + p11 * q01;
		double j10 = p00 * q10 + p10 * q11;
		double j11 = p01 * q10 + p11 * q11;
		double det = j00 * j11 - j01 * j10;
		double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;

		double det_2plus1 = 2.0 / det / det + 2.0;//4
		double det_3tr2 = 2.0*tr / det / det / det;//4
		double r0 = sD.scaf_weight * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = sD.scaf_weight * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = sD.scaf_weight * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = sD.scaf_weight * (j11 * det_2plus1 - det_3tr2 * j00);

		if (mid2varid[f0] != -1)
		{
			p_g_(scaf_F(i_p,0),0) -= (r0*p00 + r1 * p01);
			p_g_(scaf_F(i_p, 0),1) -= (r2*p00 + r3 * p01);

		}

		if (mid2varid[f1] != -1)
		{
			p_g_(scaf_F(i_p, 1),0) -= (r0*p10 + r1 * p11);
			p_g_(scaf_F(i_p, 1),1) -= (r2*p10 + r3 * p11);

		}
		if (mid2varid[f2] != -1)
		{
			p_g_(scaf_F(i_p, 2),0) += (r0*(p00 + p10) + r1 * (p01 + p11));
			p_g_(scaf_F(i_p, 2),1) += (r2*(p00 + p10) + r3 * (p01 + p11));
		}
	}

#endif // 1

	////test
	//vector<double> ppd(2 * V_N);
	//for (size_t i = 0; i < 2 * V_N; i++)
	//{
	//	ppd[i] = p_g_[i];
	//}
	//out_vec_file<double>(ppd, "d_src_grad.txt");
}

double HybridMap::max_step(const Eigen::MatrixXd& p_d_)
{
	int f0, f1, f2;
	double x02, y02, x12, y12, dx02, dy02, dx12, dy12;
	double a, b, c;
	double t_max = numeric_limits<double>::infinity();
	double t_tmp;
	block_isscaf_fid.second = -2;
	steps_.clear();
	steps_.resize(F_N + order1_fix_F_N - fix_F_N);
	for (int i = 0; i < F_N; i++)
	{
		f0 = sD.F(i, 0);
		f1 = sD.F(i, 1);
		f2 = sD.F(i, 2);
		x02 = sD.w_uv(f0,0) - sD.w_uv(f2,0);
		y02 = sD.w_uv(f0,1) - sD.w_uv(f2,1);
		x12 = sD.w_uv(f1,0) - sD.w_uv(f2,0);
		y12 = sD.w_uv(f1,1) - sD.w_uv(f2,1);
		dx02 = p_d_(f0, 0) - p_d_(f2, 0);
		dy02 = p_d_(f0, 1) - p_d_(f2, 1);
		dx12 = p_d_(f1, 0) - p_d_(f2, 0);
		dy12 = p_d_(f1, 1) - p_d_(f2, 1);

		a = dx02 * dy12 - dx12 * dy02;
		b = x02 * dy12 + dx02 * y12 - y02 * dx12 - x12 * dy02;
		c = x02 * y12 - x12 * y02;
		t_tmp = get_smallest_pos_quad_zero(a, b, c);
		steps_[i].first=t_tmp;
		steps_[i].second = i;
		if (t_max > t_tmp)
		{
			t_max = t_tmp;
			if (order1_fs.find(i) != order1_fs.end())
			{
				block_isscaf_fid.first = 0;
			}
			else
			{
				block_isscaf_fid.first = 2;
			}
			block_isscaf_fid.second = i;
		}
	}
	for (int i = fix_F_N; i < order1_fix_F_N; i++)
	{
		f0 = s_F(i, 0);
		f1 = s_F(i, 1);
		f2 = s_F(i, 2);
		x02 = s_V(f0,0) - s_V(f2,0);
		y02 = s_V(f0,1) - s_V(f2,1);
		x12 = s_V(f1,0) - s_V(f2,0);
		y12 = s_V(f1,1) - s_V(f2,1);

		int i_p = i - fix_F_N;
		f0 = scaf_F(i_p, 0);
		f1 = scaf_F(i_p, 1);
		f2 = scaf_F(i_p, 2);
		dx02 = p_d_(f0, 0) - p_d_(f2, 0);
		dy02 = p_d_(f0, 1) - p_d_(f2, 1);
		dx12 = p_d_(f1, 0) - p_d_(f2, 0);
		dy12 = p_d_(f1, 1) - p_d_(f2, 1);

		a = dx02 * dy12 - dx12 * dy02;
		b = x02 * dy12 + dx02 * y12 - y02 * dx12 - x12 * dy02;
		c = x02 * y12 - x12 * y02;
		t_tmp = get_smallest_pos_quad_zero(a, b, c);
		steps_[i].first=t_tmp;
		steps_[i_p + F_N].second = i + F_N;

		if (t_max > t_tmp)
		{
			t_max = t_tmp;
			block_isscaf_fid.first = 1;
			block_isscaf_fid.second = i;
		}
	}

	return t_max;
}

void HybridMap::backtracking_line_search(const Eigen::MatrixXd& p_d_, const Eigen::MatrixXd& neg_grad, double& step_)
{
	double h_ = 0.5;
	double c_ = 0.2;
	double tt = -(neg_grad.col(0).dot(p_d_.col(0)));
	tt-= (neg_grad.col(1).dot(p_d_.col(1)));
	//for (size_t i = 0; i < 2*pAs; i++)
	//{
	//	tt -= p_d_[i] * neg_grad[i];
	//}
	
	printf("tt==: %f\n", tt);
	Eigen::MatrixXd posAscaf;
	posAscaf.resize(pAs, 2);
	int pure_scaf_v_num_ = order1_fix_V_N - fix_V_N;

	posAscaf.topRows(V_N) = sD.w_uv;
	posAscaf.bottomRows(pure_scaf_v_num_) = s_V.block(fix_V_N, 0, pure_scaf_v_num_, 2);

	//posAscaf.block(0, 0, V_N, 1) = sD.w_uv.col(0);
	//posAscaf.block(V_N, 0, pure_scaf_v_num_, 1) = s_V.block(fix_V_N, 0, pure_scaf_v_num_, 1);
	//posAscaf.block(0, 1, V_N, 1) = sD.w_uv.col(1);
	//posAscaf.block(V_N, 1, pure_scaf_v_num_, 1) = s_V.block(fix_V_N, 1, pure_scaf_v_num_, 1);

	double ex=calc_energy(posAscaf,true);
	printf("ex==: %f\n", ex);

	Eigen::MatrixXd pos_new = posAscaf + step_ * p_d_;
	double e=calc_energy(pos_new,true);
	printf("e==: %f\n", e);

	double e_cache = e;
	double step_cache = step_;
	while (e > ex + step_ * c_*tt)
	{
		step_ *= h_;
		pos_new = posAscaf + step_ * p_d_;
		e = calc_energy(pos_new,true);
		if (e < e_cache)
		{
			e_cache = e;
			step_cache = step_;
		}
		printf("e==: %f\n", e);
	}
	//if (e_cache < e)
	//{
	//	step_ = step_cache;
	//}
}

double HybridMap::calc_energy(const Eigen::MatrixXd & pos_cur, bool with_scaf)
{
	double e_sum = 0.0;

	if (with_scaf)
	{
#pragma omp parallel for
		for (int i = 0; i < F_N; i++)
		{
			int f0 = sD.F(i, 0);
			int f1 = sD.F(i, 1);
			int f2 = sD.F(i, 2);
			double q00 = pos_cur(f0,0) - pos_cur(f2,0);
			double q01 = pos_cur(f1,0) - pos_cur(f2,0);
			double q10 = pos_cur(f0,1) - pos_cur(f2,1);
			double q11 = pos_cur(f1,1) - pos_cur(f2,1);
			//p00 = src_inv00[i]; p01 = src_inv01[i]; p10 = src_inv10[i]; p11 = src_inv11[i];
			double p00 = sD.upd_inv00[i]; 
			double p01 = sD.upd_inv01[i]; 
//			double p10 = sD.upd_inv10[i]; 
			double p11 = sD.upd_inv11[i];
			double j00 = p00 * q00;
			double j01 = p01 * q00 + p11 * q01;
			double j10 = p00 * q10;
			double j11 = p01 * q10 + p11 * q11;
			double det = j00 * j11 - j01 * j10;
			double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
			double E_d = sD.src_weight[i] * (tr + tr / (det*det));
#pragma omp critical
			{
				e_sum += E_d;
			}
		}
#pragma omp parallel for
		for (int i = 0; i < scaf_F.rows(); i++)
		{
			int f0 = scaf_F(i, 0);
			int f1 = scaf_F(i, 1);
			int f2 = scaf_F(i, 2);
			double q00 = pos_cur(f0, 0) - pos_cur(f2, 0);
			double q01 = pos_cur(f1, 0) - pos_cur(f2, 0);
			double q10 = pos_cur(f0, 1) - pos_cur(f2, 1);
			double q11 = pos_cur(f1, 1) - pos_cur(f2, 1);
			//double q00 = pos_cur(f0) - pos_cur(f2);
			//double q01 = pos_cur(f1) - pos_cur(f2);
			//double q10 = pos_cur(f0 + pAs) - pos_cur(f2 + pAs);
			//double q11 = pos_cur(f1 + pAs) - pos_cur(f2 + pAs);
			double p00 = sD.scaf_jac[i][0];
			double p01 = sD.scaf_jac[i][1];
//			double p10 = sD.scaf_jac[i][2];
			double p11 = sD.scaf_jac[i][3];
			double j00 = p00 * q00;
			double j01 = p01 * q00 + p11 * q01;
			double j10 = p00 * q10;
			double j11 = p01 * q10 + p11 * q11;
			double det = j00 * j11 - j01 * j10;
			double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
			double E_d = sD.scaf_weight * (tr + tr / (det*det));
#pragma omp critical
			{
				e_sum += E_d;
			}
		}
		return e_sum;
	}
	else
	{
		e_order_1 = pair<double, double>(0., 0.);
		e_order_non1 = pair<double, double>(0., 0.);

#pragma omp parallel for
		for (int i = 0; i < F_N; i++)
		{
			int f0 = sD.F(i, 0);
			int f1 = sD.F(i, 1);
			int f2 = sD.F(i, 2);
			double q00 = pos_cur(f0,0) - pos_cur(f2,0);
			double q01 = pos_cur(f1,0) - pos_cur(f2,0);
			double q10 = pos_cur(f0,1) - pos_cur(f2,1);
			double q11 = pos_cur(f1,1) - pos_cur(f2,1);

			double p00 = sD.upd_inv00[i];
			double p01 = sD.upd_inv01[i];
			double p10 = sD.upd_inv10[i];
			double p11 = sD.upd_inv11[i];

			//double p00 = sD.src_inv00[i]; 
			//double p01 = sD.src_inv01[i]; 
			//double p10 = sD.src_inv10[i]; 
			//double p11 = sD.src_inv11[i];
			double j00 = p00 * q00 + p10 * q01; 
			double j01 = p01 * q00 + p11 * q01; 
			double j10 = p00 * q10 + p10 * q11; 
			double j11 = p01 * q10 + p11 * q11;
			double det = j00 * j11 - j01 * j10;
			double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
			double E_d = tr + tr / (det*det);
			E_d *= sD.src_weight[i];
			sD.distortion[i] = E_d;
#pragma omp critical
			{
				bool is_boundary_region = (order1_fs.find(i) != order1_fs.end());
				if (is_boundary_region)
				{
					e_order_1.first += E_d;
					e_order_1.second += sD.src_weight[i];
				}
				else
				{
					e_order_non1.first += E_d;
					e_order_non1.second += sD.src_weight[i];
				}
				e_sum += E_d;
			}
		}


		printf("energy_order_1: %f unit %f,energy_order_non1: %f unit %f\n", e_order_1.first, e_order_1.first / e_order_1.second, e_order_non1.first, e_order_non1.first / e_order_non1.second);
		LOG(INFO) << "energy_order_1: " << e_order_1.first << " unit  " << e_order_1.first / e_order_1.second << ", energy_order_non1 : " << e_order_non1.first << " unit  " << e_order_non1.first / e_order_non1.second;
		return e_sum;
	}
}

Eigen::Vector3d HybridMap::calc_img(const Eigen::Vector3d & kesi, const std::vector<std::pair<int, BezierBase>>& bases_)
{
	double x_ = 0.0;
	double y_ = 0.0;
	for (const auto& var : bases_)
	{
		auto tri_ = var.second(kesi[0], kesi[1]);
		double coeff_ = std::get<0>(tri_);
		x_ += coeff_ * variables_[var.first];
		y_ += coeff_ * variables_[var.first + VAR_N];
	}

	return Eigen::Vector3d(x_,y_,0.0);
}

void HybridMap::engine_cond_num(set<int>& seeds_, int seeds_num)
{
//	Eigen::MatrixXd p_g_;
//	p_g_.setZero(V_N, 2);
	vector<pair<double, int>> e_dis_iso(F_N);
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
		double e_= sD.src_weight[i] * (tr + tr / (det*det));
		e_dis_iso[i].first = e_;
		e_dis_iso[i].second = i;
		sD.distortion[i] = e_;
		double det_2plus1 = 2.0 / det / det + 2.0;
		double det_3tr2 = 2.0*tr / det / det / det;
		double r0 = sD.src_weight[i] * (j00 * det_2plus1 - det_3tr2 * j11);
		double r1 = sD.src_weight[i] * (j01 * det_2plus1 + det_3tr2 * j10);
		double r2 = sD.src_weight[i] * (j10 * det_2plus1 + det_3tr2 * j01);
		double r3 = sD.src_weight[i] * (j11 * det_2plus1 - det_3tr2 * j00);

	}

	{
		set<int> fs_tmp;
		for (const auto&var : fs_tmp)
		{
			seeds_.insert(sD.F(var, 0));
		}
		printf("fs_tmp: %d ----vs_tmp: %d --\n", fs_tmp.size(),seeds_.size());
	}

	sort(e_dis_iso.begin(), e_dis_iso.end(), greater<pair<double, int>>());
	for (int i = 0; i < seeds_num; i++)
	{
		int fid_ = e_dis_iso[i].second;
		seeds_.insert(sD.F(fid_, 0));

	}
}

void HybridMap::set_minarea_for_retri(double minarea_)
{
	min_area_for_retri = minarea_;
	printf("Set the max_area for retri to %f\n", minarea_);
	LOG(INFO) << "Set the max_area for retri to %f\n" << minarea_;
}

//assembly engine_seeds_v
void HybridMap::detect_inner_engine_seeds()
{
	//add inner-v
	engine_seeds_v.clear();

	clock_t t1 = clock();
	int f_wanted_num = ((F_N / 100) > 10000 ? 10000 : F_N / 100);

	engine_cond_num(engine_seeds_v, f_wanted_num);
	clock_t t2 = clock();
	printf("detect_inner_engine_seeds----- %f -----f_wanted_num: %d ----engine_seeds_v_num: %d --\n", (t2 - t1) / 1000.0, f_wanted_num, engine_seeds_v.size());

}

void HybridMap::find_outer_rings(set<int>& part, std::vector<tuple<std::vector<int>, bool, int>>& boundaryloop)
{
	cout << "enter-------------" << endl;
	typedef tuple<int, int, int> Triplet3;
	boundaryloop.clear();
	set<Triplet3> boundary_es;

	std::vector<Triplet3> boundaryEdges;
	std::vector<Triplet3> edges;
	for (auto& var : part)
	{
		auto f_h = sD.mesh.face_handle(var);
		for (auto ifh = sD.mesh.fh_begin(f_h); ifh != sD.mesh.fh_end(f_h); ifh++)
		{
			int v0 = sD.mesh.from_vertex_handle(*ifh).idx();
			int v1 = sD.mesh.to_vertex_handle(*ifh).idx();
			if (v0 > v1) std::swap(v0, v1);
			edges.emplace_back(v0, v1, ifh->idx());
		}
	}

	std::sort(edges.begin(), edges.end());
	int i = 1;
	for (; i < edges.size();)
	{
		auto& r1 = edges[i - 1];
		auto& r2 = edges[i];
		if ((get<0>(r1) == get<0>(r2)) && (get<1>(r1) == get<1>(r2)))
		{
			i += 2;
		}
		else
		{
			boundaryEdges.emplace_back(edges[i - 1]);
			i++;
		}
	}
	if (i == edges.size())
		boundaryEdges.emplace_back(edges.back());

	for (auto&var : boundaryEdges)
	{
		auto h_h = sD.mesh.halfedge_handle(get<2>(var));
		int v0 = sD.mesh.from_vertex_handle(h_h).idx();
		int v1 = sD.mesh.to_vertex_handle(h_h).idx();
		boundary_es.emplace(v0, v1, h_h.idx());
	}

	set<int> b_es_clean;
	for (auto&var : boundary_es)
		b_es_clean.insert(get<2>(var));

	std::vector<std::vector<int>> boundaryloop_hh;
	while (!b_es_clean.empty())
	{
		int h_h_id = *(b_es_clean.begin());
		int next_h_h_id = h_h_id;
		vector<int> loop0;
		do
		{
			loop0.push_back(next_h_h_id);
			auto h_h = sD.mesh.opposite_halfedge_handle(sD.mesh.halfedge_handle(next_h_h_id));
			OpenMesh::HalfedgeHandle next_h_h = h_h;
			do
			{
				next_h_h = sD.mesh.opposite_halfedge_handle(sD.mesh.prev_halfedge_handle(next_h_h));
				if (b_es_clean.count(next_h_h.idx()) > 0)
					break;
			} while (next_h_h != h_h);
			next_h_h_id = next_h_h.idx();
			b_es_clean.erase(next_h_h_id);
		} while (next_h_h_id != h_h_id);
		boundaryloop_hh.emplace_back(loop0);
	}

	vector<int> vs_flag(V_N, -1);
	for (int i = 0; i < sD.boundary_inner_loop.size(); i++)
	{
		vs_flag[sD.boundary_inner_loop[i]] = i;
	}
	vector<tuple<int, bool, int, int,int>> b_inner_nodes;

	bool is_restart = false;
	set<int> vs_existed_;
	for (int i = 0; i < boundaryloop_hh.size(); i++)
	{
		auto& var = boundaryloop_hh[i];
		if (var.size() == 3)
		{
			int f0 = sD.mesh.face_handle(sD.mesh.opposite_halfedge_handle(sD.mesh.halfedge_handle(var[0]))).idx();
			int f1 = sD.mesh.face_handle(sD.mesh.opposite_halfedge_handle(sD.mesh.halfedge_handle(var[1]))).idx();
			int f2 = sD.mesh.face_handle(sD.mesh.opposite_halfedge_handle(sD.mesh.halfedge_handle(var[2]))).idx();
			if (f0 == f1 && f1 == f2)
			{
				part.insert(f0);
				continue;
			}
		}

		vector<int> loop0;
		int begin_v = sD.mesh.from_vertex_handle(sD.mesh.halfedge_handle(var.front())).idx();
		int end_v = sD.mesh.to_vertex_handle(sD.mesh.halfedge_handle(var.back())).idx();
		if (begin_v == end_v)
		{
			bool is_existed_ = true;
			for (auto it = var.begin(); it != var.end(); it++)
			{
				int v_id_ = sD.mesh.from_vertex_handle(sD.mesh.halfedge_handle(*it)).idx();
				loop0.push_back(v_id_);
				if (vs_existed_.find(v_id_) != vs_existed_.end())
				{
					is_restart = true;
					for (auto&itvf : sD.mesh.vf_range(sD.mesh.vertex_handle(v_id_)))
					{
						part.erase(itvf.idx());
					}
				}
				else
				{
					vs_existed_.insert(v_id_);
				}
				int fid_in_ = sD.mesh.face_handle(sD.mesh.halfedge_handle(*it)).idx();
				if (sD.order1_fs_b.find(fid_in_)== sD.order1_fs_b.end())
					is_existed_ = false;
			}
			//delete the rings existed
			if (is_existed_)
				continue;

			bool is_intersect_ = true;
			for (int j = 0; j < loop0.size(); j++)
			{
				if (vs_flag[loop0[j]] != -1)
				{
					is_intersect_=false;
					auto prev_heh_0 = sD.mesh.halfedge_handle(var[(j - 1 + var.size()) % var.size()]);
					auto prev_heh_1 = sD.mesh.opposite_halfedge_handle(prev_heh_0);
					auto next_heh_0 = sD.mesh.halfedge_handle(var[(j + var.size()) % var.size()]);
					auto next_heh_1 = sD.mesh.opposite_halfedge_handle(next_heh_0);
					bool is_out_0 = (sD.order1_fs_b.find(sD.mesh.face_handle(prev_heh_0).idx()) == sD.order1_fs_b.end()) && (sD.order1_fs_b.find(sD.mesh.face_handle(prev_heh_1).idx()) == sD.order1_fs_b.end());
					bool is_out_1 = (sD.order1_fs_b.find(sD.mesh.face_handle(next_heh_0).idx()) == sD.order1_fs_b.end()) && (sD.order1_fs_b.find(sD.mesh.face_handle(next_heh_1).idx()) == sD.order1_fs_b.end());

					if (is_out_0)
					{
						if (is_out_1)
						{
							for(auto&itvf: sD.mesh.vf_range(sD.mesh.vertex_handle(loop0[j])))
							{
								part.insert(itvf.idx());
							}
							is_restart = true;
						}
						else
						{
							//i,j,qian|false;vs_flag[loop0[j]]-1,hou|true
							b_inner_nodes.emplace_back(vs_flag[loop0[j]], true, i, j, 0);
						}
					}
					else if(is_out_1)
					{
						//i,j,hou|true;vs_flag[loop0[j]]-1,qian|false
						b_inner_nodes.emplace_back(vs_flag[loop0[j]], false, i, j, 0);
					}

					//if (fix_b_v_src2s[loop0[(j - 1 + loop0.size()) % loop0.size()]] == -1)
					//{
					//	if (fix_b_v_src2s[loop0[(j + 1) % loop0.size()]] == -1)
					//	{
					//		std::cout << "相切于Inner边界!!!" << endl;
					//	}
					//	else
					//	{
					//		//i,j,qian|false;vs_flag[loop0[j]]-1,hou|true
					//		b_inner_nodes.emplace_back(vs_flag[loop0[j]], true, i, j, 0);
					//	}
					//}
					//else if(fix_b_v_src2s[loop0[(j + 1) % loop0.size()]] == -1)
					//{
					//	//i,j,hou|true;vs_flag[loop0[j]]-1,qian|false
					//	b_inner_nodes.emplace_back(vs_flag[loop0[j]], false, i, j, 0);
					//}
				}
			}

			if (is_intersect_&&!is_restart)
			{
				bool is_counterclockwise = judge_loop_direction(loop0);
				if (is_counterclockwise)
				{
					int fid_ = sD.mesh.face_handle(sD.mesh.halfedge_handle(var.front())).idx();
					boundaryloop.emplace_back(loop0, true, fid_);
				}
				else
				{
					//在空心的地方计算空心包含的三角形个数，如果小于某个数值，则把它们全部insert到part里面。
					boundaryloop.emplace_back(loop0, false, -1);
				}
			}
		}
		else
		{
			std::cout << "ERROR===ERROR===ERROR===ERROR===ERROR=== " << var.size() << endl;
		}
	}

	if (is_restart)
	{
		boundaryloop.clear();
		find_outer_rings(part, boundaryloop);
		return;
	}

	if(b_inner_nodes.empty())
		boundaryloop.emplace_back(sD.boundary_inner_loop, false, -1);

	vector<int> nodes_flag(sD.boundary_inner_loop.size(), -1);
	for (int i = 0; i < b_inner_nodes.size(); i++)
	{
		nodes_flag[std::get<0>(b_inner_nodes[i])] = i;
	}

	for (auto&var : b_inner_nodes)
	{
		if (std::get<4>(var) == 1)
			continue;
		std::get<4>(var) = 1;
		vector<int> loop0;
		auto start_node = var;
		auto next_node = var;
		do
		{
			loop0.push_back(sD.boundary_inner_loop[std::get<0>(next_node)]);
			if (std::get<1>(next_node))
			{
				int vidx;
				for (int i = 1; i < sD.boundary_inner_loop.size(); i++)
				{
					vidx = (std::get<0>(next_node) + i) % sD.boundary_inner_loop.size();
					if (nodes_flag[vidx] != -1)
						break;
					loop0.push_back(sD.boundary_inner_loop[vidx]);
				}
				std::get<4>(b_inner_nodes[nodes_flag[vidx]]) = 1;
				next_node = b_inner_nodes[nodes_flag[vidx]];
			}
			else
			{
				vector<int> loop_;
				auto& varh = boundaryloop_hh[std::get<2>(next_node)];
				for (auto it = varh.begin(); it != varh.end(); it++)
				{
					int v_id_ = sD.mesh.from_vertex_handle(sD.mesh.halfedge_handle(*it)).idx();
					loop_.push_back(v_id_);
				}

				bool is_stayin = true;
				for (int i = 1; i < loop_.size()&&is_stayin; i++)
				{
					int vidx = (std::get<3>(next_node) + i) % loop_.size();
					for (auto&node_ : b_inner_nodes)
					{
						if (loop_[vidx] == sD.boundary_inner_loop[std::get<0>(node_)])
						{
							is_stayin = false;
							std::get<4>(node_) = 1;
							next_node = node_;
							break;
						}
					}
					if (is_stayin)
						loop0.push_back(loop_[vidx]);
				}
			}
		} while (std::get<0>(next_node) !=std::get<0>(start_node) );
		boundaryloop.emplace_back(loop0, false, -1);
	}

}

void HybridMap::set_fix_src()
{
	set<int> order1_fs_in;// = sD.patches_fids[sD.patches_cur_id];
	detect_inner_engine_seeds();
	sD.grow_from_seeds(order1_fs_in, engine_seeds_v, 5);
//	sD.find_order1_region_near_boundary_by_radius(order1_fs_in, engine_seeds_v, 5);

	//cout << "order1_fs_in: " << order1_fs_in.size() << "sD.patches_cur_id " << sD.patches_cur_id << endl;

	clock_t t1 = clock();
	boundary_loops_part.clear();
//	order1_fs_in.insert(sD.order1_fs_innerseam.begin(), sD.order1_fs_innerseam.end());
	find_outer_rings(order1_fs_in, boundary_loops_part);
	clock_t t2 = clock();
	std::cout << "find-rings" << (t2 - t1) / 1000.0 <<" order1_fs_in "<< order1_fs_in.size() << endl;

	order1_fs=sD.order1_fs_b;
	order1_fs.insert(order1_fs_in.begin(), order1_fs_in.end());
	Eigen::MatrixXi PFsrc;
	PFsrc.resize(order1_fs.size(), sD.F.cols());
	fix_f_src2s.resize(F_N);
	fix_f_src2s.setConstant(-1);
	int c_ = 0;
	for (auto&var : order1_fs)
	{
		PFsrc.row(c_) = sD.F.row(var);
		fix_f_src2s[var] = c_;
		c_++;
	}
	//fix_v_src_2_simplyMesh
	remove_unreferenced(V_N, PFsrc, fix_v_s2src, s_F, fix_v_src2s);

	fix_F_N = s_F.rows();
	fix_V_N = fix_v_s2src.size();
	printf("fix_V_N: %d , fix_F_N: %d ;\n", fix_V_N, fix_F_N);

}

bool HybridMap::judge_loop_direction(const vector<int>& loop_)
{
	double sign_area = 0.0;
	for (int i = 1; i < loop_.size(); i++)
	{
		sign_area += (sD.w_uv(loop_[i - 1],0) * sD.w_uv(loop_[i],1) - sD.w_uv(loop_[i],0) * sD.w_uv(loop_[i - 1],1));	
	}
	sign_area += (sD.w_uv(loop_.back(),0) * sD.w_uv(loop_.front(),1) - sD.w_uv(loop_.front(),0) * sD.w_uv(loop_.back(),1));

	if (sign_area > 0)
		return true;
	else
		return false;
}

void HybridMap::add_scaffold()
{
	auto& outside_b_for_retri = sD.boundary_loops.front();
	Eigen::MatrixXd bnd_pts;
	bnd_pts.resize(outside_b_for_retri.size(), 2);

	for (size_t i = 0; i < outside_b_for_retri.size(); i++)
	{
		bnd_pts.row(i) = sD.w_uv.row(outside_b_for_retri[i]);
	}
	if (sD.rect_frame_V.size() == 0) {
		Eigen::Matrix2d ob;// = rect_corners;
		{
			Eigen::VectorXd uv_max = bnd_pts.colwise().maxCoeff();
			Eigen::VectorXd uv_min = bnd_pts.colwise().minCoeff();
			Eigen::VectorXd uv_mid = (uv_max + uv_min) / 2.;


			Eigen::Array2d scaf_range(3, 3);
			ob.row(0) = uv_mid.array() + scaf_range * ((uv_min - uv_mid).array());
			ob.row(1) = uv_mid.array() + scaf_range * ((uv_max - uv_mid).array());
		}

		ob.row(1) = ob.row(1) - ob.row(0);
		for (int i = 0; i < sD.mv_num; ++i)
		{
			sD.w_uv.row(i)= sD.w_uv.row(i) - ob.row(0);
		}
		ob.row(0).setZero();

		for (size_t i = 0; i < outside_b_for_retri.size(); i++)
		{
			bnd_pts.row(i) = sD.w_uv.row(outside_b_for_retri[i]);
		}

		Eigen::Vector2d rect_len;
		rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);
		int frame_points = 20;
		sD.interval = rect_len(0);

		sD.rect_frame_V.resize(4 * frame_points, 2);
		for (int i = 0; i < frame_points; i++) {
			// 0,0;0,1
			sD.rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
			// 0,0;1,1
			sD.rect_frame_V.row(i + frame_points)
				<< ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
			// 1,0;1,1
			sD.rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
				- i * rect_len(1) /
				frame_points;
			// 1,0;0,1
			sD.rect_frame_V.row(i + 3 * frame_points)
				<< ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
			// 0,0;0,1
		}	
		sD.frame_ids = Eigen::VectorXi::LinSpaced(sD.rect_frame_V.rows(), fix_V_N,
			fix_V_N + sD.rect_frame_V.rows());
	}

	Eigen::MatrixXd retri_s_V;
	Eigen::MatrixXi retri_s_E;

	retri_s_V.resize(bnd_pts.rows() + sD.rect_frame_V.rows(), 2);
	retri_s_V << bnd_pts, sD.rect_frame_V;

	retri_s_E.resize(retri_s_V.rows(), 2);
	for (int i = 0; i < retri_s_E.rows(); i++)
		retri_s_E.row(i) << i, i + 1;
	int bnd_size_ = bnd_pts.rows();
	retri_s_E(bnd_size_ - 1, 1) = 0;
	retri_s_E(retri_s_V.rows() - 1, 1) = bnd_size_;

	Eigen::MatrixXd retri_s_H;
	retri_s_H.setZero(sD.component_sizes.size(), 2);
	{
		int hole_f = 0;
		int hole_i = 0;
		for (auto cs : sD.component_sizes) 
		{
			for (int i = 0; i < 3; i++)
			{
				retri_s_H.row(hole_i) += sD.w_uv.row(sD.F(hole_f, i));
			}
			hole_f += cs;
			hole_i++;
		}
	}
	retri_s_H /= 3.;

	Eigen::MatrixXd retri_s_uv;
	triangulate(retri_s_V, retri_s_E, retri_s_H, retri_s_uv, sD.s_T);
	cout << "add scaffold: vnums " << retri_s_uv.rows() << " fnums " << sD.s_T.rows() << endl;

	for (int i = 0; i < sD.s_T.rows(); i++)
	{
		for (int j = 0; j < sD.s_T.cols(); j++)
		{
			auto& v_id = sD.s_T(i, j);
			if (v_id < outside_b_for_retri.size())
				v_id = fix_v_src2s[outside_b_for_retri[v_id]];
			else
				v_id += (fix_V_N - outside_b_for_retri.size());
		}
	}

	s_F.conservativeResize(fix_F_N + sD.s_T.rows(), 3);
	s_F.bottomRows(sD.s_T.rows()) = sD.s_T;
	s_V.resize(fix_V_N + retri_s_uv.rows()- bnd_pts.rows(), 2);
	for (size_t i = 0; i < fix_V_N; i++)
	{
		s_V.row(i) = sD.w_uv.row(fix_v_s2src[i]);
	}
	s_V.bottomRows(retri_s_uv.rows() - bnd_pts.rows()) = retri_s_uv.bottomRows(retri_s_uv.rows() - bnd_pts.rows());

	order1_fix_F_N = s_F.rows();
	order1_fix_V_N = s_V.rows();
	printf("order1_fix_V_N: %d , order1_fix_F_N: %d ;\n", order1_fix_V_N, order1_fix_F_N);

//	writeObj(s_V, s_F, "ij.obj");
}

void HybridMap::reTriangulation()
{
	printf("Boundary_loops_part has %d loops.\n", boundary_loops_part.size());
	if (boundary_loops_part.empty())
		return;
	Eigen::MatrixXd retri_in_H;
	Eigen::MatrixXi retri_in_E;
	Eigen::MatrixXd pts;
	Eigen::MatrixXi T_F;
	Eigen::MatrixXd bnd_pts;

	scafd.internal_bnd.resize(0);
	scafd.bnd_sizes.clear();
	for (auto&var : boundary_loops_part)
	{
		auto& loop_ = std::get<0>(var);
		scafd.internal_bnd.conservativeResize(scafd.internal_bnd.size() + loop_.size());
		scafd.internal_bnd.bottomRows(loop_.size()) = Eigen::Map<Eigen::ArrayXi>(loop_.data(), loop_.size());
		scafd.bnd_sizes.push_back(loop_.size());
		if (std::get<1>(var))
		{
			retri_in_H.conservativeResize(retri_in_H.rows() + 1, 2);
			int fid_ = std::get<2>(var);
			//Eigen::Vector2d hole_=(sD.w_uv.row(sD.F(fid_, 0)) + sD.w_uv.row(sD.F(fid_, 1)) + sD.w_uv.row(sD.F(fid_, 2))) / 3.0;
			retri_in_H.bottomRows(1) = (sD.w_uv.row(sD.F(fid_, 0)) + sD.w_uv.row(sD.F(fid_, 1)) + sD.w_uv.row(sD.F(fid_, 2))) / 3.0;
		}		
	}

	bnd_pts.resize(scafd.internal_bnd.size(), 2);	
	for (size_t i = 0; i < scafd.internal_bnd.size(); i++)
	{
		bnd_pts.row(i) = sD.w_uv.row(scafd.internal_bnd[i]);
	}

	retri_in_E.resize(bnd_pts.rows(), 2);
	for (int i = 0; i < bnd_pts.rows(); i++)
	{
		retri_in_E.row(i) << i, i + 1;
	}
	int c_ = 0;
	for (const auto&var : scafd.bnd_sizes)
	{
		retri_in_E(c_ + var - 1, 1) = c_;
		c_ += var;
	}

	triangulate(bnd_pts, retri_in_E, retri_in_H, min_area_for_retri, pts, T_F);

	for (size_t i = 0; i < T_F.rows(); i++)
	{
		for (size_t j = 0; j < T_F.cols(); j++)
		{
			auto& v_id = T_F(i, j);
			if (v_id < scafd.internal_bnd.size())
				v_id = fix_v_src2s[scafd.internal_bnd[v_id]];
			else
			{
				v_id += (order1_fix_V_N - bnd_pts.rows());
			}
		}
	}
	int inner_num_without_segLine = pts.rows() - bnd_pts.rows();
	s_V.conservativeResize(order1_fix_V_N + inner_num_without_segLine, 2);
	s_V.bottomRows(inner_num_without_segLine) = pts.bottomRows(inner_num_without_segLine);

	s_F.conservativeResize(order1_fix_F_N + T_F.rows(), 3);
	s_F.bottomRows(T_F.rows()) = T_F;

	//pos_of_Smesh.resize(s_V.rows() * 2);
	printf("After reTriangulation, there are %d faces and %d vertices;\n", s_F.rows(), s_V.rows());
	LOG(INFO) << "fix_F_N: " << fix_F_N << " order1_fix_F_N: " << order1_fix_F_N << " total F_N: " << s_F.rows();
	LOG(INFO) << "fix_V_N: " << fix_V_N << " order1_fix_V_N: " << order1_fix_V_N << " total V_N: " << s_V.rows();

//	writeObj(s_V, s_F, "ij2.obj");

}

void HybridMap::meshid_vs_varid()
{
	varid2mid.clear();
	varid2mid.reserve(s_V.rows());
	for (int i = 0; i < fix_V_N; i++)
	{
		varid2mid.push_back(i);
	}
	for (int i = fix_V_N+sD.frame_ids.size(); i < s_V.rows(); i++)
	{
		varid2mid.push_back(i);
	}

	mid2varid.assign(s_V.rows(), -1);
	for (int i = 0; i < varid2mid.size(); i++)
	{
		mid2varid[varid2mid[i]] = i;
	}
}

//prepare V_N,F_N,V,F,VF;allocate memory for pos_of_mesh, distortion;
//rescale the mesh surface area to pi
void HybridMap::prepare_data()
{	
	V_N = sD.mv_num;
	F_N = sD.mf_num;
	fine2coarse_v_id_pos.resize(V_N);

	samples_fid_sim.resize(F_N);
	samples_bc.resize(F_N);

//	pos_of_mesh.resize(2 * V_N);

	if (!is_init_2d)
		update_pos_from_sD();

	block_isscaf_fid.first = -1;

}

//update pos_of_mesh, update mesh for show; calc the boundary_loops
void HybridMap::tutte_map()
{
	Eigen::MatrixXd uv_init;
	for (auto&var : sD.boundary_loops)
	{
		Eigen::MatrixXd bnd_uv;
		Eigen::VectorXi bnd = Eigen::Map<Eigen::VectorXi>(var.data(), var.size());
		map_vertices_to_circle(sD.V, bnd, bnd_uv);
		bnd_uv *= sqrt(1.0 / M_PI);
		Tutte(V_N, sD.F, bnd, bnd_uv, uv_init);
	}

	sD.w_uv = uv_init;
}

void HybridMap::run_1_iter(int iter_num)
{
	clock_t t1, t2, t_begin, t_end ,tti;
	vector<double> t_vecc;
	t_vecc.reserve(20);
	t2 = clock();
	t_begin = clock();
	t_vecc.emplace_back(clock() / 1000.0);
	set_fix_src();
	t_vecc.emplace_back(clock() / 1000.0);

	add_scaffold();
	t_vecc.emplace_back(clock() / 1000.0);

	min_area_for_retri = 0.01;


	//cycle operator
	reTriangulation();
	t_vecc.emplace_back(clock() / 1000.0);


	meshid_vs_varid();
	t_vecc.emplace_back(clock() / 1000.0);

	//return;
	gen_s_mesh(s_V, s_F);
	t_vecc.emplace_back(clock() / 1000.0);

	t1 = clock();

	//printf("gen-s-mesh finish----- %f ------sth begin-----------\n", (t1 - t2) / 1000.0);
	find_index_fine2coarse();
	t_vecc.emplace_back(clock() / 1000.0);

	t2 = clock();


	calc_scaf_jacobian();
	t_vecc.emplace_back(clock() / 1000.0);

	setOrder();
	t_vecc.emplace_back(clock() / 1000.0);

	t1 = clock();
	//printf("calc-jac && set-order finish----- %f ------sth begin-----------\n", (t1 - t2) / 1000.0);

	setNodes();
	setIdentity();

	t_vecc.emplace_back(clock() / 1000.0);

	t_end = clock();
	t2 = clock();

	solve_equation();
	t_vecc.emplace_back(clock() / 1000.0);

	printf("set_fix_src//add_scaffold//reTriangulation//meshid_vs_varid//gen_s_mesh//find_index_fine2coarse//calc_scaf_jacobian//setOrder//setNodes&setIdentity//solve_equation \n");
	for (int ti = 0; ti < t_vecc.size() - 1; ti++)
	{
		cout << t_vecc[ti + 1] - t_vecc[ti] << "//";
	}
	cout << endl;

}

void HybridMap::autorun(double conv_rate_)
{ 
	enter_this_solver();
	int iter_num_cur = 0;
	double energy_pre;
	double conv_percent=1.0;
	while (iter_num_cur < MAX_ITER_NUM)
	{
		printf("Iteration begin-------------------------%d------------------------\n",iter_num_cur);
		LOG(INFO) << "Iteration begin-------------------------"<< iter_num_cur <<"------------------------";
		energy_pre = sD.Energy_cur;
		if(sD.Energy_cur>20)
			is_slim_meth = true;
		else
		{
			is_slim_meth = false;
		}
		clock_t t_b = clock();
		run_1_iter(iter_num_cur);
		clock_t t_e = clock();
		printf("Iteration end--------------------%d-----------time %f-------------\n", iter_num_cur,(t_e-t_b)/1000.0);
		LOG(INFO) << "Iteration begin-------------------------" << iter_num_cur << "------------time "<< (t_e - t_b) / 1000.0 <<"-------------";
		sD.energy_record.push_back(sD.Energy_cur);
		sD.energy_record3.emplace_back(sD.Energy_cur, (t_e - t_b) / 1000.0, 0);
		iter_num_cur++;

		conv_percent = abs(energy_pre - sD.Energy_cur) / energy_pre;

		if (sD.Energy_cur > 10)
		{
			if (conv_percent <= convergence_rate)
				break;
		}else if (conv_percent <= conv_rate_ ||(sD.Energy_cur<4.8&&conv_percent<1e-2))
		{
			break;
		}
	}
	leave_this_solver();

	//local_opt_with_fixed_boundary();
}

void HybridMap::setOrder()
{
	//tri_has_samps.clear();
	//tri_has_samps.resize(s_F.rows());
	vector<int>	flag_tri_num(s_F.rows(), 0);
	vector<double> flag_tri_dis(s_F.rows(), 0.0);
	for (int i = 0; i < F_N; i++)
	{
		if (fix_f_src2s[i] != -1)
		{
			flag_tri_num[samples_fid_sim[i].front()] = sample_num_per_tri;
		}
		else
		{
			for (int j = 0; j < sample_num_per_tri; j++)
			{
				flag_tri_num[samples_fid_sim[i][j]] += 1;
				flag_tri_dis[samples_fid_sim[i][j]] += sD.distortion[i];
				//tri_has_samps[fine2coarse_f[i][j]].emplace_back(i, j);
			}
		}
	}

	//out_vec_file<int>(flag_tri_num, "tri_num.txt");
	//out_vec_file<double>(flag_tri_dis, "tri_dis.txt");

	//[ 0 =1= 100 =2= 1000 =3= ]
	flag_order.clear();
	flag_order.resize(s_F.rows(), 1);
	vector<int> ods_dis(3, 0);


	if (order_set1 == 0)
	{
		for (size_t i = 0; i < s_F.rows(); i++)
		{
			if (flag_tri_num[i] < 100)
			{
				flag_order[i] = 1;
				ods_dis[0]++;
			}
			else if (flag_tri_num[i] < 1000)
			{
				flag_order[i] = 2;
				ods_dis[1]++;
			}
			else
			{
				flag_order[i] = 3;
				ods_dis[2]++;
			}
		}
	}
	else if(order_set1 == 3)
	{
		flag_order.clear();
		flag_order.resize(s_F.rows(), 3);
		for (int i = 0; i < order1_fix_F_N; i++)
		{
			flag_order[i] = 1;
		}
		ods_dis[0] = order1_fix_F_N;
		ods_dis[2] = s_F.rows() - order1_fix_F_N;
	}

	std::printf("order1 = %d; order2 = %d; order3 = %d;\n", ods_dis[0], ods_dis[1], ods_dis[2]);
	LOG(INFO) << "order1 = " << ods_dis[0] << "; order2 = " << ods_dis[1] << "; order3 = " << ods_dis[2];
}

//find the location index of the barycenter of fineMesh triangle, store them to fine2coarse_f 
void HybridMap::find_index_fine2coarse()
{
	//calc area of s_mesh
	vector<double> area_s_mesh(s_F.rows());
#pragma omp parallel for
	for (int i = 0; i < s_F.rows(); i++)
	{
		int f0 = s_F(i, 0);
		int f1 = s_F(i, 1);
		int f2 = s_F(i, 2);
		double p00 = s_V(f1, 0) - s_V(f0, 0);
		double p10 = s_V(f1, 1) - s_V(f0, 1);
		double p01 = s_V(f2, 0) - s_V(f0, 0);
		double p11 = s_V(f2, 1) - s_V(f0, 1);
		area_s_mesh[i] = p00 * p11 - p01 * p10;
	}
	FindIndex findidx(s_V, s_F);

	//calc for samples
#pragma omp parallel for
	for (int i = 0; i < F_N; i++)
	{
		//cout << i << endl;
		if (fix_f_src2s[i] != -1)
		{
			samples_fid_sim[i] = vector<int>{ fix_f_src2s[i] };
			samples_bc[i] = vector<Eigen::Vector2d>{ Eigen::Vector2d(1. / 3., 1. / 3.) };
		}
		else
		{
			vector<int> f_sample_(sample_num_per_tri);
			vector<Eigen::Vector2d> f_s_pos;
			f_s_pos.reserve(sample_num_per_tri);
			Eigen::Vector3d bc_;
			int idx_;
			Eigen::Vector3d p0(sD.w_uv(sD.F(i, 0),0), sD.w_uv(sD.F(i, 0),1), 0.0);
			Eigen::Vector3d p1(sD.w_uv(sD.F(i, 1),0), sD.w_uv(sD.F(i, 1),1), 0.0);
			Eigen::Vector3d p2(sD.w_uv(sD.F(i, 2),0), sD.w_uv(sD.F(i, 2),1), 0.0);

			bc_ = (p0 + p1 + p2) / 3.0;
			idx_ = findidx.query_point(bc_);
			f_sample_[0] = idx_;
			calc_bc(bc_, idx_, area_s_mesh);
			f_s_pos.emplace_back(bc_[0],bc_[1]);
			/////////////////////////////////
			bc_ = (p0 + 4.0*p1 + p2) / 6.0;
			idx_ = findidx.query_point(bc_);
			f_sample_[1] = idx_;
			calc_bc(bc_, idx_, area_s_mesh);
			f_s_pos.emplace_back(bc_[0], bc_[1]);
			/////////////////////////////////
			bc_ = (p0 + p1 + 4.0*p2) / 6.0;
			idx_ = findidx.query_point(bc_);
			f_sample_[2] = idx_;
			calc_bc(bc_, idx_, area_s_mesh);
			f_s_pos.emplace_back(bc_[0], bc_[1]);
			/////////////////////////////////
			bc_ = (4.0*p0 + p1 + p2) / 6.0;
			idx_ = findidx.query_point(bc_);
			f_sample_[3] = idx_;
			calc_bc(bc_, idx_, area_s_mesh);
			f_s_pos.emplace_back(bc_[0], bc_[1]);
			/////////////////////////////////

			samples_fid_sim[i] = f_sample_;
			samples_bc[i] = f_s_pos;
		}
	}

	//calc for source vertex 
#pragma omp parallel for
	for (int i = 0; i < V_N; i++)
	{
		if (fix_v_src2s[i] != -1)
			continue;
		Eigen::Vector3d bc_(sD.w_uv(i,0), sD.w_uv(i,1), 0.0);
		int idx_ = findidx.query_point(bc_);
		calc_bc(bc_, idx_, area_s_mesh);
	//	fine2coarse_vid_sim[i] = idx_;
	//	fine2coarse_v_bc[i] = Eigen::Vector2d(bc_[0], bc_[1]);

		fine2coarse_v_id_pos[i] = make_pair(idx_, Eigen::Vector2d(bc_[0], bc_[1]));

	}

}

void HybridMap::gen_s_mesh(const Eigen::MatrixXd & V_, const Eigen::MatrixXi & F_)
{
	s_mesh.clear();
	int dim_v = V_.cols();
	if (dim_v == 2)
	{
		for (size_t i = 0; i < V_.rows(); i++)
		{
			OpenMesh::Vec3d pos_(V_(i, 0), V_(i, 1), 0.0);
			s_mesh.add_vertex(pos_);
		}
	}
	else if (dim_v == 3)
	{
		for (size_t i = 0; i < V_.rows(); i++)
		{
			OpenMesh::Vec3d pos_(V_(i, 0), V_(i, 1), V_(i,2));
			s_mesh.add_vertex(pos_);
		}
	}
	for (size_t i = 0; i < F_.rows(); i++)
	{
		std::vector<OpenMesh::VertexHandle> face_v_hs;
		for (size_t j = 0; j < F_.cols(); j++)
		{
			face_v_hs.push_back(s_mesh.vertex_handle(F_(i, j)));
		}
		s_mesh.add_face(face_v_hs);
	}
	printf("simply-mesh-----vertex: %d -----faces: %d -----\n", s_mesh.n_vertices(),s_mesh.n_faces());
	//OpenMesh::IO::write_mesh(s_mesh, "jj.obj");
	//writeObj(s_V, s_F, "zijixie.obj");
}

void HybridMap::setNodes()
{
	//variable V

	int var_v_num = s_V.rows() - sD.frame_ids.size();
	int count_ = var_v_num;
	for (auto ite = s_mesh.edges_begin(); ite != s_mesh.edges_end(); ite++)
	{
		auto he_h = s_mesh.halfedge_handle(*ite, 0);
		auto he_h_op = s_mesh.halfedge_handle(*ite, 1);
		auto f0_h = s_mesh.face_handle(he_h);
		auto f1_h = s_mesh.face_handle(s_mesh.opposite_halfedge_handle(he_h));
		if (flag_order[f0_h.idx()] == 1 || flag_order[f1_h.idx()] == 1)
		{
			continue;
		}
		else
		{
			int order_sum = flag_order[f0_h.idx()] + flag_order[f1_h.idx()];
			if (order_sum <= 5)
			{
				s_mesh.data(he_h).var_ids.push_back(count_);
				s_mesh.data(he_h_op).var_ids.push_back(count_);
				count_++;
			}
			else
			{
				s_mesh.data(he_h).var_ids.push_back(count_);
				s_mesh.data(he_h_op).var_ids.push_back(count_ + 1);
				count_++;
				s_mesh.data(he_h).var_ids.push_back(count_);
				s_mesh.data(he_h_op).var_ids.push_back(count_ - 1);
				count_++;
			}
		}
	}
	int var_e_num = count_ - var_v_num;

	for (auto itf = s_mesh.faces_begin(); itf != s_mesh.faces_end(); itf++)
	{
		auto he_h2 = *(s_mesh.fh_begin(*itf));
		int v_h0_idx = mid2varid[s_mesh.to_vertex_handle(he_h2).idx()];
		auto he_h0 = s_mesh.next_halfedge_handle(he_h2);
		int v_h1_idx = mid2varid[s_mesh.to_vertex_handle(he_h0).idx()];
		auto he_h1 = s_mesh.next_halfedge_handle(he_h0);
		int v_h2_idx = mid2varid[s_mesh.to_vertex_handle(he_h1).idx()];
		auto& f_ids_bases = s_mesh.data(*itf).var_ids_bases;

		if (itf->idx() < fix_F_N)
		{
			f_ids_bases.emplace_back(v_h0_idx, bezier_100);
			f_ids_bases.emplace_back(v_h1_idx, bezier_010);
			f_ids_bases.emplace_back(v_h2_idx, bezier_001);
		}
		else if (itf->idx() < order1_fix_F_N)
		{
			if (v_h0_idx != -1)
				f_ids_bases.emplace_back(v_h0_idx, bezier_100);
			if (v_h1_idx != -1)
				f_ids_bases.emplace_back(v_h1_idx, bezier_010);
			if (v_h2_idx != -1)
				f_ids_bases.emplace_back(v_h2_idx, bezier_001);
		}
		else
		{
			if (flag_order[itf->idx()] == 1)
			{
				f_ids_bases.emplace_back(v_h0_idx, bezier_100);
				f_ids_bases.emplace_back(v_h1_idx, bezier_010);
				f_ids_bases.emplace_back(v_h2_idx, bezier_001);
			}
			else if (flag_order[itf->idx()] == 2)
			{
				if (s_mesh.data(he_h2).var_ids.size() == 1)
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_101);
					if (s_mesh.data(he_h0).var_ids.size() == 1)
					{
						f_ids_bases.emplace_back(v_h0_idx, bezier_200);
						f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_110);
						if (s_mesh.data(he_h1).var_ids.size() == 1)
						{
							//222
							f_ids_bases.emplace_back(v_h1_idx, bezier_020);
							f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_011);
							f_ids_bases.emplace_back(v_h2_idx, bezier_002);
						}
						else
						{
							//221
							f_ids_bases.emplace_back(v_h1_idx, bezier_2_211);
							f_ids_bases.emplace_back(v_h2_idx, bezier_2_122);
						}
					}
					else
					{
						f_ids_bases.emplace_back(v_h0_idx, bezier_2_201);
						if (s_mesh.data(he_h1).var_ids.size() == 1)
						{
							//212
							f_ids_bases.emplace_back(v_h1_idx, bezier_2_112);
							f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_011);
							f_ids_bases.emplace_back(v_h2_idx, bezier_002);
						}
						else
						{
							//211
							f_ids_bases.emplace_back(v_h1_idx, bezier_2_111);
							f_ids_bases.emplace_back(v_h2_idx, bezier_2_122);
						}
					}
				}
				else
				{
					if (s_mesh.data(he_h0).var_ids.size() == 1)
					{
						f_ids_bases.emplace_back(v_h0_idx, bezier_2_102);
						f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_110);
						if (s_mesh.data(he_h1).var_ids.size() == 1)
						{
							//122
							f_ids_bases.emplace_back(v_h1_idx, bezier_020);
							f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_011);
							f_ids_bases.emplace_back(v_h2_idx, bezier_2_221);
						}
						else
						{
							//121
							f_ids_bases.emplace_back(v_h1_idx, bezier_2_211);
							f_ids_bases.emplace_back(v_h2_idx, bezier_2_121);
						}
					}
					else
					{
						f_ids_bases.emplace_back(v_h0_idx, bezier_2_101);
						if (s_mesh.data(he_h1).var_ids.size() == 1)
						{
							//112
							f_ids_bases.emplace_back(v_h1_idx, bezier_2_112);
							f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_011);
							f_ids_bases.emplace_back(v_h2_idx, bezier_2_221);
						}
						else
						{
							//111
							f_ids_bases.emplace_back(v_h1_idx, bezier_2_111);
							f_ids_bases.emplace_back(v_h2_idx, bezier_2_121);
						}
					}
				}
			}
			else if (flag_order[itf->idx()] == 3)
			{
				int i_e2 = s_mesh.data(he_h2).var_ids.size();
				int i_e0 = s_mesh.data(he_h0).var_ids.size();
				int i_e1 = s_mesh.data(he_h1).var_ids.size();

				int type_i = i_e2 * 100 + i_e0 * 10 + i_e1;
				switch (type_i)
				{
				case 222://333
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_300);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_030);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_003);
				}
				break;
				case 221://332
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_300);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_312);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_223);
				}
				break;
				case 220://331
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_300);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_311);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_123);
				}
				break;
				case 212://323
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_302);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_213);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_003);
				}
				break;
				case 211://322
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_302);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_212);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_223);
				}
				break;
				case 210://321
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_302);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_211);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_123);
				}
				break;
				case 202://313
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_301);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_113);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_003);
				}
				break;
				case 201://312
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_301);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_112);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_223);
				}
				break;
				case 200://311
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.back(), bezier_201);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_301);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_111);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_123);
				}
				break;

				//
				case 122://233
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_203);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_030);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_322);
				}
				break;
				case 121://232
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_203);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_312);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_222);
				}
				break;
				case 120://231
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_203);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_311);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_122);
				}
				break;
				case 112://223
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_202);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_213);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_322);
				}
				break;
				case 111://222
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_202);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_212);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_222);
				}
				break;
				case 110://221
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_202);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_211);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_122);
				}
				break;
				case 102://213
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_201);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_113);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_322);
				}
				break;
				case 101://212
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_201);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_112);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_222);
				}
				break;
				case 100://211
				{
					f_ids_bases.emplace_back(s_mesh.data(he_h2).var_ids.front(), bezier_3_2_e20);
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_201);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_111);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_122);
				}
				break;

				//
				case 22://133
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_103);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_030);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_321);
				}
				break;
				case 21://132
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_103);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_312);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_221);
				}
				break;
				case 20://131
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_103);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_210);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.back(), bezier_120);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_311);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_121);
				}
				break;
				case 12://123
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_213);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_321);
				}
				break;
				case 11://122
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_212);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_221);
				}
				break;
				case 10://121
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_102);
					f_ids_bases.emplace_back(s_mesh.data(he_h0).var_ids.front(), bezier_3_2_e01);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_211);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_121);
				}
				break;
				case 2://113
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_101);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_113);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_021);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.back(), bezier_012);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_321);
				}
				break;
				case 1://112
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_101);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_112);
					f_ids_bases.emplace_back(s_mesh.data(he_h1).var_ids.front(), bezier_3_2_e12);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_221);
				}
				break;
				case 0://111
				{
					f_ids_bases.emplace_back(v_h0_idx, bezier_3_101);
					f_ids_bases.emplace_back(v_h1_idx, bezier_3_111);
					f_ids_bases.emplace_back(v_h2_idx, bezier_3_121);
				}
				break;
				default:
					break;
				}
				s_mesh.data(*itf).var_ids_bases.emplace_back(count_, bezier_111);
				count_++;
			}
		}
	}

	int var_f_num = count_ - var_v_num - var_e_num;
	VAR_N = count_;
	printf("var_num: %d === var_v_num: %d; var_e_num: %d; var_f_num: %d;\n", count_, var_v_num, var_e_num, var_f_num);
	LOG(INFO) << "var_num: "<< count_ <<" === var_v_num: "<< var_v_num <<"; var_e_num: "<< var_e_num <<"; var_f_num: "<< var_f_num <<";";
}

void HybridMap::setIdentity()
{
	variables_.resize(2 * VAR_N);

	for (int itv = 0; itv < s_V.rows(); itv++)
	{
		if (mid2varid[itv] != -1)
		{
			variables_[mid2varid[itv]] = s_V(itv, 0);
			variables_[mid2varid[itv] + VAR_N] = s_V(itv, 1);
		}

	}
	for (auto ite = s_mesh.edges_begin(); ite != s_mesh.edges_end(); ite++)
	{
		auto he_h = s_mesh.halfedge_handle(*ite, 0);
		auto he_h_op = s_mesh.halfedge_handle(*ite, 1);
		auto f0_h = s_mesh.face_handle(he_h);
		auto f1_h = s_mesh.face_handle(s_mesh.opposite_halfedge_handle(he_h));

		auto vf_h = s_mesh.from_vertex_handle(he_h);
		auto vt_h = s_mesh.to_vertex_handle(he_h);
		auto vf_pos = s_mesh.point(vf_h);
		auto vt_pos = s_mesh.point(vt_h);

		if (flag_order[f0_h.idx()] == 1 || flag_order[f1_h.idx()] == 1)
		{
			continue;
		}
		else
		{
			int order_sum = flag_order[f0_h.idx()] + flag_order[f1_h.idx()];
			if (order_sum <= 5)
			{
				int cur_id=s_mesh.data(he_h).var_ids.front();
				OpenMesh::Vec3d vm_pos = (vf_pos + vt_pos) / 2.0;
				variables_[cur_id] = vm_pos[0];
				variables_[cur_id + VAR_N] = vm_pos[1];
			}
			else
			{
				int cur_id = s_mesh.data(he_h).var_ids.front();
				OpenMesh::Vec3d vm_pos = (2.0*vf_pos + vt_pos) / 3.0;
				variables_[cur_id] = vm_pos[0];
				variables_[cur_id + VAR_N] = vm_pos[1];
				cur_id = s_mesh.data(he_h).var_ids.back();
				vm_pos = (vf_pos + 2.0*vt_pos) / 3.0;
				variables_[cur_id] = vm_pos[0];
				variables_[cur_id + VAR_N] = vm_pos[1];
			}
		}
	}

	for (auto itf = s_mesh.faces_begin(); itf != s_mesh.faces_end(); itf++)
	{
		if (flag_order[itf->idx()] == 3)
		{
			int cur_id = s_mesh.data(*itf).var_ids_bases.back().first;
			OpenMesh::Vec3d vpos_(0., 0., 0.);
			for (auto&fv_h : s_mesh.fv_range(*itf))
			{
				vpos_ += s_mesh.point(fv_h);
			}
			vpos_ /= 3.0;
			variables_[cur_id] = vpos_[0];
			variables_[cur_id + VAR_N] = vpos_[1];
		}
	}


}

void HybridMap::assembly_CM_para()
{
	vector<double>& gradient_ = pardiso_solver->rhs;
	vector<double>& ps_a = pardiso_solver->a;
	ps_a.clear();
	gradient_.clear();
	ps_a.resize(pardiso_solver->ja.size(), 0.0);
	gradient_.resize(VAR_N * 2, 0.0);

	//vector<double> e_b_tmp(F_N*sample_num_per_tri, 0.0);

#pragma omp parallel for
	for (int trii = 0; trii < F_N; trii++)
	{
		//cout << trii << endl;
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		hess_compents_record.reserve(sample_num_per_tri * 300);
		grad_compents_record.reserve(sample_num_per_tri * 50);

		for (int trij = 0; trij < sample_num_per_tri; trij++)
		{
			int i = samples_fid_sim[trii][trij];
			auto& fh = s_mesh.face_handle(i);
			auto& ids_bases = s_mesh.data(fh).var_ids_bases;
			int bases_num = ids_bases.size();
			bool is_boundary_region = (i < order1_fix_F_N);

			double f_area = (is_boundary_region ? sD.src_weight[trii] : sD.src_weight[trii] / sample_num_per_tri);

			double p00, p01, p10, p11;
			if (is_boundary_region)
			{
				p00 = sD.upd_inv00[trii];
				p01 = sD.upd_inv01[trii];
				p10 = sD.upd_inv10[trii];
				p11 = sD.upd_inv11[trii];
			}
			else
			{
				int f0 = sD.F(trii, 0);
				int f1 = sD.F(trii, 1);
				int f2 = sD.F(trii, 2);
				double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
				double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
				double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
				double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);
				double j00 = q00 * sD.upd_inv00[trii] + q01 * sD.upd_inv10[trii];
				double j01 = q00 * sD.upd_inv01[trii] + q01 * sD.upd_inv11[trii];
				double j10 = q10 * sD.upd_inv00[trii] + q11 * sD.upd_inv10[trii];
				double j11 = q10 * sD.upd_inv01[trii] + q11 * sD.upd_inv11[trii];
				int f0_s = s_F(i, 0);
				int f1_s = s_F(i, 1);
				int f2_s = s_F(i, 2);
				double s00 = s_V(f1_s, 1) - s_V(f2_s, 1);
				double s01 = s_V(f2_s, 0) - s_V(f1_s, 0);
				double s10 = s_V(f2_s, 1) - s_V(f0_s, 1);
				double s11 = s_V(f0_s, 0) - s_V(f2_s, 0);
				double det_ = s00 * s11 - s01 * s10;
				p00 = (s00 * j00 + s01 * j10) / det_;
				p01 = (s00 * j01 + s01 * j11) / det_;
				p10 = (s10 * j00 + s11 * j10) / det_;
				p11 = (s10 * j01 + s11 * j11) / det_;
			}


			const auto& bc_ = samples_bc[trii][trij];
			vector<double> der_kesi0(bases_num);
			vector<double> der_kesi1(bases_num);
			for (size_t j = 0; j < bases_num; j++)
			{
				auto& res_ = ids_bases[j].second(bc_[0], bc_[1]);
				der_kesi0[j] = std::get<1>(res_);
				der_kesi1[j] = std::get<2>(res_);
			}
			vector<double> coeff_ac(bases_num);
			vector<double> coeff_bd(bases_num);
			double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
			for (size_t j = 0; j < bases_num; j++)
			{
				coeff_ac[j] = p00 * der_kesi0[j] + p10 * der_kesi1[j];
				coeff_bd[j] = p01 * der_kesi0[j] + p11 * der_kesi1[j];
				a += variables_[ids_bases[j].first] * coeff_ac[j];
				b += variables_[ids_bases[j].first] * coeff_bd[j];
				c += variables_[ids_bases[j].first + VAR_N] * coeff_ac[j];
				d += variables_[ids_bases[j].first + VAR_N] * coeff_bd[j];
			}


			double tr = a * a + b * b + c * c + d * d;
			double det = a * d - b * c;
			double alpha = sqrt(tr + 2 * det);
			double beta = sqrt((a - d)*(a - d) + (b + c)*(b + c));
			double sig1 = (alpha + beta) * 0.5;
			double sig2 = (alpha - beta) * 0.5;
			double det_1 = 1.0 / det;
			double det_2plus1 = det_1 * det_1;
			double tr_det_4_3 = 3.*tr*det_2plus1 * det_2plus1;
			double det_3_2 = 2.*det_2plus1 * det_1;

			double det_3tr2 = det_2plus1 * tr*det_1;
			det_2plus1 += 1.0;
			double r0 = f_area * (a * det_2plus1 - det_3tr2 * d);
			double r1 = f_area * (b * det_2plus1 + det_3tr2 * c);
			double r2 = f_area * (c * det_2plus1 + det_3tr2 * b);
			double r3 = f_area * (d * det_2plus1 - det_3tr2 * a);

			//calc gradient
			for (size_t j = 0; j < bases_num; j++)
			{
				double g_u_ = r0 * coeff_ac[j] + r1 * coeff_bd[j];
				double g_v_ = r2 * coeff_ac[j] + r3 * coeff_bd[j];
				grad_compents_record.emplace_back(ids_bases[j].first, g_u_);
				grad_compents_record.emplace_back(ids_bases[j].first + VAR_N, g_v_);
			}


			double uu0 = (det_2plus1 + tr_det_4_3 *d*d - 2.*det_3_2*a*d)*f_area;
			double uu12 = (det_3_2*(a*c - b * d) - tr_det_4_3*d*c)*f_area;
			double uu3 = (det_2plus1 + tr_det_4_3*c*c + 2.*det_3_2*b*c)*f_area;

			double vv0 = (det_2plus1 + tr_det_4_3*b*b + 2.*det_3_2*b*c)*f_area;
			double vv12 = (det_3_2*(b * d- a * c) - tr_det_4_3*a*b)*f_area;
			double vv3 = (det_2plus1 + tr_det_4_3*a*a - 2.*det_3_2*a*d)*f_area;

			double uv0 = (det_3_2*(a*b - c * d) - tr_det_4_3*b*d)*f_area;
			double uv1 = (tr_det_4_3*a*d - det_3_2*(a*a + d * d + tr / 2.))*f_area;
			double uv2 = (tr_det_4_3*b*c + det_3_2*(b*b + c * c + tr / 2.))*f_area;
			double uv3 = (det_3_2*(c * d- a * b) - tr_det_4_3*a*c)*f_area;

		//	double sig_flag = (sig1 + sig2 - 1.0 / sig1 / sig1 / sig1 - 1.0 / sig2 / sig2 / sig2);
			double sig_flag = (sig1 + sig2 - 1.0 / (sig1 * sig1 * sig1) - 1.0 / (sig2 * sig2 * sig2));

			double sig_u0, sig_u1, sig_u2;
			bool is_sig_flag = (sig_flag < 0);
			if(is_sig_flag)
			{
				double alpha_1 = 1.0 / alpha;
				double alpha_3 = 0.5*alpha_1 * alpha_1*alpha_1;
				sig_u0 = (b - c)*(b - c)*alpha_3;
				sig_u1 = (a + d)*(b - c)*alpha_3;
				sig_u2 = (a + d)*(a + d)*alpha_3;
				sig_flag *= f_area;
				sig_u0 *= sig_flag;
				sig_u1 *= sig_flag;
				sig_u2 *= sig_flag;

				//sig_u0//-sig_u1//..//sig_u2
				//sig_u2//sig_u1//..//sig_u0

				//sig_u1//sig_u0//-sig_u2//-sig_u1

			}

			//calc hessian
			for (size_t ui = 0; ui < bases_num; ui++)
			{
				for (size_t uj = 0; uj < bases_num; uj++)
				{

					double uij11 = coeff_ac[ui] * coeff_ac[uj];
					double uij12 = coeff_ac[ui] * coeff_bd[uj];
					double uij21 = coeff_bd[ui] * coeff_ac[uj];
					double uij22 = coeff_bd[ui] * coeff_bd[uj];

					double h_ui_vj_plus = uv0 * uij11 + uv1 * uij12 + uv2 * uij21 + uv3 * uij22;
					if (is_sig_flag)
					{
						h_ui_vj_plus -= (sig_u1 * uij11 + sig_u0 * uij12 - sig_u2 * uij21 - sig_u1 * uij22);
					}
					hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first] + find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first + VAR_N), h_ui_vj_plus);
					if (ids_bases[uj].first >= ids_bases[ui].first)
					{
						double h_ui_uj_plus = uu0 * uij11 + uu12 * uij12 + uu12 * uij21 + uu3 * uij22;
						double h_vi_vj_plus = vv0 * uij11 + vv12 * uij12 + vv12 * uij21 + vv3 * uij22;
						if (is_sig_flag)
						{
							h_ui_uj_plus -= (sig_u0 * uij11 - sig_u1 * uij12 - sig_u1 * uij21 + sig_u2 * uij22);
							h_vi_vj_plus -= (sig_u2 * uij11 + sig_u1 * uij12 + sig_u1 * uij21 + sig_u0 * uij22);
						}
						int dd_ = find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first);
						hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first] + dd_, h_ui_uj_plus);
						hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first + VAR_N] + dd_, h_vi_vj_plus);
					}
				}
			}

			if (is_boundary_region)
				break;
		}

#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				gradient_[var.first] -= var.second;
		}
	}

	//for (int i = 0; i < core_num; i++)
	//{
	//	for (const auto&var : hess_compents_record[i])
	//		ps_a[var.first] += var.second;
	//	for (const auto&varg : grad_compents_record[i])
	//		gradient_[varg.first] -= varg.second;
	//}

	double scaf_w_ = sD.Energy_cur / (order1_fix_F_N - fix_F_N) / 100.0;
	printf("energy_beginning: %f,scaffold-weight: %f\n", sD.Energy_cur, scaf_w_);
	sD.adjust_scaf_weight(scaf_w_);
	assembly_scaf_part_para();

	////test
	//vector<double> ppd(2 * V_N);
	//for (int i = 0; i < V_N; i++)
	//{
	//	int s_vid = fix_v_src2s[i];
	//	if (s_vid != -1)
	//	{
	//		ppd[i] = gradient_[s_vid];
	//		ppd[i + V_N] = gradient_[s_vid + VAR_N];
	//	}
	//	else
	//	{
	//	}
	//}
	//out_vec_file<double>(ppd, "grad.txt");

}

void HybridMap::assembly_SLIM_para()
{
	vector<double>& gradient_ = pardiso_solver->rhs;
	vector<double>& ps_a = pardiso_solver->a;
	ps_a.clear();
	gradient_.clear();
	ps_a.resize(pardiso_solver->ja.size(), 0.0);
	gradient_.resize(VAR_N * 2, 0.0);

	//double weight_adjust = e_order_non1.second / e_order_1.second;

#pragma omp parallel for
	for (int trii = 0; trii < F_N; trii++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		hess_compents_record.reserve(sample_num_per_tri * 300);
		grad_compents_record.reserve(sample_num_per_tri * 50);

		for (int trij = 0; trij < sample_num_per_tri; trij++)
		{
			int i = samples_fid_sim[trii][trij];
			auto& fh = s_mesh.face_handle(i);
			auto& ids_bases = s_mesh.data(fh).var_ids_bases;
			int bases_num = ids_bases.size();
			bool is_boundary_region = (i < order1_fix_F_N);
			
			double f_area = (is_boundary_region ? sD.src_weight[trii] : sD.src_weight[trii] / sample_num_per_tri);

			double p00, p01, p10, p11;
			if (is_boundary_region)
			{
				p00 = sD.upd_inv00[trii];
				p01 = sD.upd_inv01[trii];
				p10 = sD.upd_inv10[trii];
				p11 = sD.upd_inv11[trii];
			}
			else
			{
				int f0 = sD.F(trii, 0);
				int f1 = sD.F(trii, 1);
				int f2 = sD.F(trii, 2);
				double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
				double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
				double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
				double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);
				double j00 = q00 * sD.upd_inv00[trii] + q01 * sD.upd_inv10[trii];
				double j01 = q00 * sD.upd_inv01[trii] + q01 * sD.upd_inv11[trii];
				double j10 = q10 * sD.upd_inv00[trii] + q11 * sD.upd_inv10[trii];
				double j11 = q10 * sD.upd_inv01[trii] + q11 * sD.upd_inv11[trii];
				int f0_s = s_F(i, 0);
				int f1_s = s_F(i, 1);
				int f2_s = s_F(i, 2);
				double s00 = s_V(f1_s, 1) - s_V(f2_s, 1);
				double s01 = s_V(f2_s, 0) - s_V(f1_s, 0);
				double s10 = s_V(f2_s, 1) - s_V(f0_s, 1);
				double s11 = s_V(f0_s, 0) - s_V(f2_s, 0);
				double det_ = s00 * s11 - s01 * s10;
				p00 = (s00 * j00 + s01 * j10)/det_;
				p01 = (s00 * j01 + s01 * j11)/det_;
				p10 = (s10 * j00 + s11 * j10)/det_;
				p11 = (s10 * j01 + s11 * j11)/det_;
			}

			//const auto& c_jac = const_jac[trii][trij];
			//double p00 = c_jac[0];
			//double p01 = c_jac[1];
			//double p10 = c_jac[2];
			//double p11 = c_jac[3];

			const auto& bc_ = samples_bc[trii][trij];
			vector<double> der_kesi0(bases_num);
			vector<double> der_kesi1(bases_num);
			for (size_t j = 0; j < bases_num; j++)
			{
				auto& res_ = ids_bases[j].second(bc_[0], bc_[1]);
				der_kesi0[j] = std::get<1>(res_);
				der_kesi1[j] = std::get<2>(res_);
			}
			vector<double> coeff_ac(bases_num);
			vector<double> coeff_bd(bases_num);
			double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
			for (size_t j = 0; j < bases_num; j++)
			{
				coeff_ac[j] = p00 * der_kesi0[j] + p10 * der_kesi1[j];
				coeff_bd[j] = p01 * der_kesi0[j] + p11 * der_kesi1[j];
				a += variables_[ids_bases[j].first] * coeff_ac[j];
				b += variables_[ids_bases[j].first] * coeff_bd[j];
				c += variables_[ids_bases[j].first + VAR_N] * coeff_ac[j];
				d += variables_[ids_bases[j].first + VAR_N] * coeff_bd[j];
			}

			/*a_u=coeff_ac, b_u=coeff_bd, a_v=b_v=0,
			  c_v=coeff_ac, d_v=coeff_bd, c_u=d_u=0,
			  tr_u,tr_v
			  det_u,det_v
			  alpha_u,alpha_v
			*/

			double tr = a * a + b * b + c * c + d * d;
			double det = a * d - b * c;
			double alpha = sqrt(tr + 2 * det);
			double beta = sqrt((a - d)*(a - d) + (b + c)*(b + c));
			double sig1 = (alpha + beta) / 2.0;
			double sig2 = (alpha - beta) / 2.0;
			double det_1 = 1.0 / det;
			double det_2plus1 = det_1*det_1;
			double det_3tr2 = det_2plus1 * tr*det_1;
			det_2plus1 += 1.0;
			double r0 = f_area * (a * det_2plus1 - det_3tr2 * d);
			double r1 = f_area * (b * det_2plus1 + det_3tr2 * c);
			double r2 = f_area * (c * det_2plus1 + det_3tr2 * b);
			double r3 = f_area * (d * det_2plus1 - det_3tr2 * a);


			double sig1_1 = 1.0 / sig1;
			double sig2_1 = 1.0 / sig2;

			double sig1_new = sqrt(1.0 + sig1_1 + sig1_1 * sig1_1 + sig1_1 * sig1_1*sig1_1);
			double sig2_new = sqrt(1.0 + sig2_1 + sig2_1 * sig2_1 + sig2_1 * sig2_1*sig2_1);

			double tmp_ = 0;
			if (beta > 1e-10)
				tmp_ = (sig2_new - sig1_new) / (sig2*sig2 - sig1 * sig1);

			double sig_square_sum = 0.5*(sig2 * sig2 + sig1 * sig1);
			double sig_new_sum = 0.5*(sig2_new + sig1_new);

			double w00 = tmp_ * (a*a + b * b - sig_square_sum) + sig_new_sum;
			double w01 = tmp_ * (a*c + b * d);
			double w11 = tmp_ * (c*c + d * d - sig_square_sum) + sig_new_sum;

			double w1 = f_area * (w00 * w00 + w01 * w01);
			double w2 = f_area * (w00 + w11)* w01;
			double w3 = f_area * (w01 * w01 + w11 * w11);

			//calc gradient
			for (size_t j = 0; j < bases_num; j++)
			{
				double g_u_ = r0*coeff_ac[j]+r1*coeff_bd[j];
				double g_v_ = r2*coeff_ac[j]+r3*coeff_bd[j];
				grad_compents_record.emplace_back(ids_bases[j].first, g_u_);
				grad_compents_record.emplace_back(ids_bases[j].first + VAR_N, g_v_);
			}

			//calc hessian
			for (size_t ui = 0; ui < bases_num; ui++)
			{
				for (size_t uj = 0; uj < bases_num; uj++)
				{
					double h_ui_uj_plus = coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj];

					hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first]+ find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first + VAR_N), w2*h_ui_uj_plus);

					if (ids_bases[uj].first >= ids_bases[ui].first)
					{
						int dd_ = find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first);
						hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first] + dd_, w1*h_ui_uj_plus);
						hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first+VAR_N] + dd_, w3*h_ui_uj_plus);
					}
				}
			}

			if (is_boundary_region)
				break;
		}

#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				gradient_[var.first] -= var.second;
		}
	}

	

	//double energy_beginning = 0.0;
	//for (const double& var : e_b_tmp)
	//	energy_beginning += var;

	double scaf_w_ = sD.Energy_cur / (order1_fix_F_N - fix_F_N) / 100.0;
	printf("energy_beginning: %f,scaffold-weight: %f\n", sD.Energy_cur, scaf_w_);

	//double scaf_w_ = (weight_adjust*e_order_1.first+e_order_non1.first) / (order1_fix_F_N - fix_F_N) / 100.0;
	//printf("energy_beginning:order1 %f, order-non1: %f, sum: %f, scaffold-weight: %f\n", weight_adjust*e_order_1.first, e_order_non1.first, (weight_adjust*e_order_1.first + e_order_non1.first),scaf_w_);

	sD.adjust_scaf_weight(scaf_w_);
	assembly_scaf_part_SLIM_para();

}

void HybridMap::assembly_scaf_part_SLIM_para()
{
	vector<double>& ps_a = pardiso_solver->a;
//	vector <double>& gradient_ = pardiso_solver->rhs;
#pragma omp parallel for
	for (int trii = fix_F_N; trii < order1_fix_F_N; trii++)
	{
		vector<pair<int, double>> hess_compents_record;
		hess_compents_record.reserve(21);
		int i = trii - fix_F_N;
		auto& fh = s_mesh.face_handle(trii);
		auto& ids_bases = s_mesh.data(fh).var_ids_bases;
		int bases_num = ids_bases.size();
		if (bases_num == 0)
			continue;

		double f_area = sD.scaf_weight;
		double p00 = sD.scaf_jac[i][0];
		double p01 = sD.scaf_jac[i][1];
		double p10 = sD.scaf_jac[i][2];
		double p11 = sD.scaf_jac[i][3];
		vector<double> der_kesi0(bases_num);
		vector<double> der_kesi1(bases_num);
		for (size_t j = 0; j < bases_num; j++)
		{
			auto& res_ = ids_bases[j].second(0.33, 0.33);
			der_kesi0[j] = std::get<1>(res_);
			der_kesi1[j] = std::get<2>(res_);
		}
		vector<double> coeff_ac(bases_num);
		vector<double> coeff_bd(bases_num);
		for (size_t j = 0; j < bases_num; j++)
		{
			coeff_ac[j] = p00 * der_kesi0[j] + p10 * der_kesi1[j];
			coeff_bd[j] = p01 * der_kesi0[j] + p11 * der_kesi1[j];
		}

		double w13 = f_area * 4.;
		//calc hessian
		for (size_t ui = 0; ui < bases_num; ui++)
		{
			for (size_t uj = 0; uj < bases_num; uj++)
			{
				if (ids_bases[uj].first >= ids_bases[ui].first)
				{
					double h_ui_uj_plus = coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj];
					int dd_ = find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first);
					hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first] + dd_, w13*h_ui_uj_plus);
					hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first + VAR_N] + dd_, w13*h_ui_uj_plus);

				}
			}
		}
#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_a[var.first] += var.second;
		}
	}
}

void HybridMap::assembly_scaf_part_para()
{
	vector<double>& ps_a = pardiso_solver->a;
	vector <double>& gradient_ = pardiso_solver->rhs;

#pragma omp parallel for
	for (int trii = fix_F_N; trii < order1_fix_F_N; trii++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		hess_compents_record.reserve(21);
		grad_compents_record.reserve(6);

		int i = trii-fix_F_N;
		auto& fh = s_mesh.face_handle(trii);
		auto& ids_bases = s_mesh.data(fh).var_ids_bases;
		int bases_num = ids_bases.size();
		if (bases_num == 0)
			continue;

		double f_area = sD.scaf_weight;
		double p00 = sD.scaf_jac[i][0];
		double p01 = sD.scaf_jac[i][1];
		double p10 = sD.scaf_jac[i][2];
		double p11 = sD.scaf_jac[i][3];
		vector<double> der_kesi0(bases_num);
		vector<double> der_kesi1(bases_num);
		for (size_t j = 0; j < bases_num; j++)
		{
			auto& res_ = ids_bases[j].second(0.33, 0.33);
			der_kesi0[j] = std::get<1>(res_);
			der_kesi1[j] = std::get<2>(res_);
		}
		vector<double> coeff_ac(bases_num);
		vector<double> coeff_bd(bases_num);
		for (size_t j = 0; j < bases_num; j++)
		{
			coeff_ac[j] = p00 * der_kesi0[j] + p10 * der_kesi1[j];
			coeff_bd[j] = p01 * der_kesi0[j] + p11 * der_kesi1[j];
		}
		int f0 = s_F(trii, 0);
		int f1 = s_F(trii, 1);
		int f2 = s_F(trii, 2);
		double q00 = s_V(f0, 0) - s_V(f2, 0);
		double q01 = s_V(f1, 0) - s_V(f2, 0);
		double a = p00 * q00 + p10 * q01;
		double b = p01 * q00 + p11 * q01;

		//a=d   b=-c
		f_area *= 2.;

		double a2p1 = (1. + a * a)*f_area;
		double ab = a * b*f_area;
		double b2p1 = (1. + b * b)*f_area;
		double uv1 = (a*a - 1.)*f_area;
		double uv2 = (1. - b * b)*f_area;
		// a2p1 ab b2p1
		// b2p1 -ab a2p1
		//-ab uv1 uv2 ab

		//double uu0 = (1. + a * a)*f_area;
		//double uu12 = (a*b)*f_area;
		//double uu3 = (1. + b*b)*f_area;

		//double vv0 = (1. + b*b)*f_area;
		//double vv12 = (- a*b)*f_area;
		//double vv3 = (1. + a*a)*f_area;

		//double uv0 = (- b*a)*f_area;
		//double uv1 = (a*a - 1.)*f_area;
		//double uv2 = (1. - b * b)*f_area;
		//double uv3 = (a*b)*f_area;

		//calc hessian
		for (size_t ui = 0; ui < bases_num; ui++)
		{
			for (size_t uj = 0; uj < bases_num; uj++)
			{

				double uij11 = coeff_ac[ui] * coeff_ac[uj];
				double uij12 = coeff_ac[ui] * coeff_bd[uj];
				double uij21 = coeff_bd[ui] * coeff_ac[uj];
				double uij22 = coeff_bd[ui] * coeff_bd[uj];

				double h_ui_vj_plus = -ab * uij11 + uv1 * uij12 + uv2 * uij21 + ab * uij22;

				hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first] + find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first + VAR_N), h_ui_vj_plus);

				if (ids_bases[uj].first >= ids_bases[ui].first)
				{
					double h_ui_uj_plus = a2p1 * uij11 + ab * uij12 + ab * uij21 + b2p1 * uij22;
					double h_vi_vj_plus = b2p1 * uij11 - ab * uij12 - ab * uij21 + a2p1 * uij22;

					int dd_ = find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first);
					hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first] + dd_, h_ui_uj_plus);
					hess_compents_record.emplace_back(pardiso_solver->ia[ids_bases[ui].first + VAR_N] + dd_, h_vi_vj_plus);
				}
			}
		}


#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				gradient_[var.first] -= var.second;
		}
	}

}

void HybridMap::calc_scaf_jacobian()
{
	//calc the area_threshold;
	double min_bnd_edge_len = numeric_limits<double>::infinity();

	auto& outside_b_for_retri = sD.boundary_loops.front();
	for (size_t i = 0; i < outside_b_for_retri.size()-1; i++)
	{
		min_bnd_edge_len = min(min_bnd_edge_len, (s_V.row(fix_v_src2s[outside_b_for_retri[i]]) - s_V.row(fix_v_src2s[outside_b_for_retri[i + 1]])).squaredNorm());
	}
	min_bnd_edge_len = min(min_bnd_edge_len, (s_V.row(fix_v_src2s[outside_b_for_retri.front()]) - s_V.row(fix_v_src2s[outside_b_for_retri.back()])).squaredNorm());
	double area_threshold = min_bnd_edge_len/4.0;
	sD.scaf_jac.resize(order1_fix_F_N - fix_F_N);
#pragma omp parallel for
	for (int i = fix_F_N; i < order1_fix_F_N; i++)
	{
		int f0 = s_F(i, 0);
		int f1 = s_F(i, 1);
		int f2 = s_F(i, 2);

		Eigen::Vector3d x_(s_V(f0, 0) - s_V(f2, 0), s_V(f0, 1) - s_V(f2, 1), 0.);
		Eigen::Vector3d l_(s_V(f1, 0) - s_V(f2, 0), s_V(f1, 1) - s_V(f2, 1), 0.);
		Eigen::Vector3d n_ = x_.cross(l_);
		double area_f = n_.norm() / 2.0;

		double x1_0, x2_0, y2_0;

		if (area_f > area_threshold)
		{
			x1_0 = x_.norm();
			x_ /= x1_0;
			n_ /= (2.0*area_f);
			Eigen::Vector3d y_ = n_.cross(x_);
			x2_0 = l_.dot(x_);
			y2_0 = l_.dot(y_);
		}
		else
		{
			cout << "find wrong!!!!!!!!!!!!!!!!!!!!!!" << endl;
			double h = sqrt((2 * area_threshold) / sqrt(3.0));
			x1_0 = h;
			x2_0 = h / 2.0;
			y2_0 = sqrt(3.0)*h / 2.0;
		}
		sD.scaf_jac[i - fix_F_N] = Eigen::Vector4d(1.0 / x1_0, -x2_0 / (x1_0*y2_0), 0.0, 1.0 / y2_0);
	}

}

void HybridMap::prepare_pardiso()
{
	if (pardiso_solver != NULL)
	{
		delete pardiso_solver;
		pardiso_solver = NULL;
	}
	pardiso_solver = new PardisoSolver();

	auto& ps_ia = pardiso_solver->ia;
	auto& ps_ja = pardiso_solver->ja;

	vector<set<int>> vv_infl(VAR_N * 2);
	for (auto itf = s_mesh.faces_begin(); itf != s_mesh.faces_end(); itf++)
	{
		auto& ids_bases = s_mesh.data(*itf).var_ids_bases;
		int bases_num = ids_bases.size();
		for (size_t i = 0; i < bases_num; i++)
		{
			for (size_t j = 0; j < bases_num; j++)
			{
				vv_infl[ids_bases[i].first].insert(ids_bases[j].first);
			}
		}
	}

	ps_ia.clear(); ps_ia.reserve(2 * VAR_N + 1);
	ps_ja.clear(); ps_ja.reserve(36 * VAR_N);

	typedef Eigen::Triplet<int> T;
	std::vector<T> tripletlist;
	tripletlist.reserve(72 * VAR_N);
	for (size_t i = 0; i < VAR_N; i++)
	{
		ps_ia.push_back(ps_ja.size());
		vector<int> row_id(vv_infl[i].begin(), vv_infl[i].end());
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			ps_ja.push_back(row_id[k]);
			tripletlist.push_back(T(i, row_id[k], dd));
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			ps_ja.push_back(row_id[k] + VAR_N);
			tripletlist.push_back(T(i, row_id[k] + VAR_N, dd));
			++dd;
		}
	}
	for (size_t i = VAR_N; i < 2*VAR_N; i++)
	{
		ps_ia.push_back(ps_ja.size());
		vector<int> row_id(vv_infl[i - VAR_N].begin(), vv_infl[i - VAR_N].end());
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i-VAR_N);

		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			ps_ja.push_back(row_id[k] + VAR_N);
		}
	}
	ps_ia.push_back(ps_ja.size());

	find_idx_in_hess.resize(VAR_N, 2 * VAR_N);
	find_idx_in_hess.setFromTriplets(tripletlist.begin(), tripletlist.end());
	//find_idx_in_hess.setFromTriplets(tripletlist.begin(), tripletlist.end());
	
	pardiso_solver->nnz = ps_ja.size();
	pardiso_solver->num = 2 * VAR_N;
	cout << "actually hessian non-zero-num: " << pardiso_solver->nnz << " actually variable V_N: " << VAR_N << endl;
	LOG(INFO) << "actually hessian non-zero-num: " << pardiso_solver->nnz << " actually variable V_N: " << VAR_N;

}

void HybridMap::calc_bc(Eigen::Vector3d& bc_, const int& f_id, const vector<double>& areas)
{
	Eigen::Vector2d p0 = s_V.row(s_F(f_id, 0));
	Eigen::Vector2d p1 = s_V.row(s_F(f_id, 1));
	Eigen::Vector2d p2 = s_V.row(s_F(f_id, 2));

	double p00 = p1(0)-p0(0);
	double p10 = p1(1)-p0(1);
	double p01 = bc_(0)-p0(0);
	double p11 = bc_(1)-p0(1);
	bc_[2] = (p00 * p11 - p01 * p10) / areas[f_id];
	//
	p00 = p2(0) - p1(0);
	p10 = p2(1) - p1(1);
	p01 = bc_(0) - p1(0);
	p11 = bc_(1) - p1(1);
	bc_[0] = (p00 * p11 - p01 * p10) / areas[f_id];
	bc_[1] = 1.0 - bc_[0] - bc_[2];

}


void HybridMap::update_pos_from_sD()
{
	sD.w_uv.conservativeResize(V_N, 2);
}

void HybridMap::engine_convergence_points(const vector<pair<double, int>>& e_dis_iso, set<int>& seeds_, int loops_num)
{
	printf("engine_convergence_points begin--------------\n");
	set<int> fs_set;
	vector<int> is_excluded(F_N, 0);

#pragma omp parallel for
	for (int i = 0; i < F_N; i++)
	{
		if (is_excluded[i] == 1)
			continue;
		vector<int> f_visit(F_N, 0);
		vector<int> v_visit(V_N, 0);
		queue<pair<int, int>> Q;
		Q.emplace(i, 0);
		f_visit[i] = 1;
		bool stop_ = false;
		while (!Q.empty()&&!stop_)
		{
			auto seed_ = Q.front();
			Q.pop();
			if (seed_.second >= loops_num)
				continue;
			for (size_t j = 0; j < sD.F.cols(); j++)
			{
				int fvid_ = sD.F(seed_.first, j);
				if (v_visit[fvid_] == 0)
				{
					v_visit[fvid_] = 1;
					auto v_h = sD.mesh.vertex_handle(fvid_);
					for (auto itvf = sD.mesh.vf_begin(v_h); itvf != sD.mesh.vf_end(v_h); itvf++)
					{
						int fid_ = itvf->idx();
						if (f_visit[fid_] == 0)
						{
							f_visit[fid_] = 1;
							if (e_dis_iso[fid_].first < e_dis_iso[i].first)
							{
								Q.emplace(fid_, seed_.second + 1);
								is_excluded[fid_] = 1;
							}
							else
							{
								is_excluded[i] = 1;
								stop_ = true;
								break;
							}
						}
					}
					//for (const auto& fid_ : sD.VF[fvid_])
					//{
					//	if (f_visit[fid_] == 0)
					//	{
					//		f_visit[fid_] = 1;
					//		if (e_dis_iso[fid_].first < e_dis_iso[i].first)
					//		{
					//			Q.emplace(fid_, seed_.second + 1);
					//			is_excluded[fid_] = 1;
					//		}
					//		else
					//		{
					//			is_excluded[i] = 1;
					//			stop_ = true;
					//			break;
					//		}
					//	}
					//}
					if (stop_)
						break;
				}
			}
		}
		if (!stop_)
		{
			is_excluded[i] = 0;
			fs_set.insert(i);
		}
	}

	int ex_count = 0;
	for (size_t i = 0; i < F_N; i++)
	{
		if (is_excluded[i] == 0)
		{
			cout << i << " ";
			ex_count++;
		}
	}
	cout << endl << "========================================================="<<ex_count << endl;
	for (auto&var : fs_set)
		cout << var << " ";
	cout << endl << "=========================================================" <<fs_set.size()<< endl;

	for (auto&fid_ : fs_set)
	{
		seeds_.insert(sD.F(fid_, 0));
		seeds_.insert(sD.F(fid_, 1));
		seeds_.insert(sD.F(fid_, 2));
	}

}

void HybridMap::load_2d_mesh(const char * filename_)
{
	Mesh mesh_2d_;
	OpenMesh::IO::read_mesh(mesh_2d_, filename_);
	int v_2d_N = mesh_2d_.n_vertices();
	int f_2d_N = mesh_2d_.n_faces();
	sD.w_uv.resize(v_2d_N, 2);
	for (auto itv = mesh_2d_.vertices_begin(); itv != mesh_2d_.vertices_end(); itv++)
	{
		auto pos_ = mesh_2d_.point(*itv);
		sD.w_uv(itv->idx(),0) = pos_[0];
		sD.w_uv(itv->idx(),1) = pos_[1];
	}
	is_init_2d = true;
	printf("Load-2d-Mesh.......................V: %d;F: %d;\n", v_2d_N, f_2d_N);
}

void HybridMap::leave_this_solver()
{
	if (pardiso_solver != NULL)
	{
		delete pardiso_solver;
		pardiso_solver = NULL;
	}
	find_idx_in_hess.resize(0, 0);
	//fine2coarse_vid_sim.clear();
	//fine2coarse_v_bc.clear();
	fine2coarse_v_id_pos.clear();
	samples_bc.clear();
	samples_fid_sim.clear();
	boundary_loops_part.clear();

	fix_v_src2s.resize(0);
	fix_v_s2src.resize(0);
	fix_f_src2s.resize(0);

	s_mesh.clear();
	s_V.resize(0,0);
	s_F.resize(0,0);
	flag_order.clear();

	variables_.resize(0);
	scaf_F.resize(0,0);
	engine_seeds_v.clear();
}

void HybridMap::enter_this_solver()
{
	//fine2coarse_vid_sim.resize(F_N);
	//fine2coarse_v_bc.resize(F_N);
	fine2coarse_v_id_pos.resize(F_N);
	samples_fid_sim.resize(F_N);
	samples_bc.resize(F_N);

	//pos_of_mesh.resize(2 * V_N);
}

void HybridMap::info_src_hessian()
{
	int src_hess_nnz=preCalc_pardiso(sD.V, sD.F);
	LOG(INFO) << "src hessian non-zero-num: " << src_hess_nnz << " src variable V_N: " << sD.mv_num;
}

pair<int, int> HybridMap::info_sim_hessian()
{
	enter_this_solver();

//	is_slim_meth = true;

	set_fix_src();

	add_scaffold();

	//cycle operator
	reTriangulation();


	meshid_vs_varid();

	//return;
	gen_s_mesh(s_V, s_F);

	find_index_fine2coarse();


	calc_scaf_jacobian();

	setOrder();

	setNodes();
	setIdentity();


	vector<int> ps_ia;
	vector<int> ps_ja;

	vector<set<int>> vv_infl(VAR_N * 2);
	for (auto itf = s_mesh.faces_begin(); itf != s_mesh.faces_end(); itf++)
	{
		auto& ids_bases = s_mesh.data(*itf).var_ids_bases;
		int bases_num = ids_bases.size();
		for (size_t i = 0; i < bases_num; i++)
		{
			for (size_t j = 0; j < bases_num; j++)
			{
				vv_infl[ids_bases[i].first].insert(ids_bases[j].first);
			}
		}
	}

	ps_ia.clear(); ps_ia.reserve(2 * VAR_N + 1);
	ps_ja.clear(); ps_ja.reserve(36 * VAR_N);

	typedef Eigen::Triplet<int> T;
	std::vector<T> tripletlist;
	tripletlist.reserve(72 * VAR_N);
	for (size_t i = 0; i < VAR_N; i++)
	{
		ps_ia.push_back(ps_ja.size());
		vector<int> row_id(vv_infl[i].begin(), vv_infl[i].end());
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			ps_ja.push_back(row_id[k]);
			tripletlist.push_back(T(i, row_id[k], dd));
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			ps_ja.push_back(row_id[k] + VAR_N);
			tripletlist.push_back(T(i, row_id[k] + VAR_N, dd));
			++dd;
		}
	}
	for (size_t i = VAR_N; i < 2 * VAR_N; i++)
	{
		ps_ia.push_back(ps_ja.size());
		vector<int> row_id(vv_infl[i - VAR_N].begin(), vv_infl[i - VAR_N].end());
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i - VAR_N);

		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			ps_ja.push_back(row_id[k] + VAR_N);
		}
	}
	ps_ia.push_back(ps_ja.size());

	int nnz_ = ps_ja.size();
	cout << "actually hessian non-zero-num: " << nnz_ << " actually variable V_N: " << VAR_N << endl;
	LOG(INFO) << "actually hessian non-zero-num: " << nnz_ << " actually variable V_N: " << VAR_N;
	int var_n_=VAR_N;

	leave_this_solver();

	return pair<int, int>(nnz_, var_n_);
}
