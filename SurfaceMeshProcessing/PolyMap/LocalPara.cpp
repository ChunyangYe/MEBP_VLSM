#include "LocalPara.h"


LocalPara::LocalPara(SharedData & shardd_, const set<int>& fids_p_, vector<int>& f_visit_):sD(shardd_),f_visit(f_visit_)
{
	v_src_N = sD.V_num;
	fids_p.assign(fids_p_.begin(), fids_p_.end());
}

LocalPara::~LocalPara()
{
}

void LocalPara::assembly_pardiso_with_fixed_boundary()
{
	ps_local.a.clear();
	ps_local.a.assign(ps_local.ja.size(), 0.0);
	ps_local.rhs.clear();
	ps_local.rhs.assign(2 * v_free_N, 0.0);

#pragma omp parallel for
	for (int trii = 0; trii < F_p.rows(); trii++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		hess_compents_record.reserve(21);
		grad_compents_record.reserve(6);

		int f0 = F_p(trii, 0);
		int f1 = F_p(trii, 1);
		int f2 = F_p(trii, 2);

		double f_area = weight_p[trii];
		double p00 = src_inv_p[trii][0];
		double p01 = src_inv_p[trii][1];
		double p10 = src_inv_p[trii][2];
		double p11 = src_inv_p[trii][3];
		vector<double> coeff_ac;
		vector<double> coeff_bd;
		vector<int>ids_bases;
		if (f0 < v_free_N)
		{
			coeff_ac.push_back(p00);
			coeff_bd.push_back(p01);
			ids_bases.push_back(f0);
		}
		if (f1 < v_free_N)
		{
			coeff_ac.push_back(p10);
			coeff_bd.push_back(p11);
			ids_bases.push_back(f1);
		}
		if (f2 < v_free_N)
		{
			coeff_ac.push_back(-(p00 + p10));
			coeff_bd.push_back(-(p01 + p11));
			ids_bases.push_back(f2);
		}
		int bases_num = ids_bases.size();

		if (bases_num == 0)
			continue;

		double q00 = V_p(f0,0) - V_p(f2,0);
		double q10 = V_p(f0,1) - V_p(f2,1);
		double q01 = V_p(f1,0) - V_p(f2,0);
		double q11 = V_p(f1,1) - V_p(f2,1);
		double a = p00 * q00 + p10 * q01;
		double b = p01 * q00 + p11 * q01;
		double c = p00 * q10 + p10 * q11;
		double d = p01 * q10 + p11 * q11;

		/*a_u=coeff_ac, b_u=coeff_bd, a_v=b_v=0,
		  c_v=coeff_ac, d_v=coeff_bd, c_u=d_u=0,
		  tr_u,tr_v
		  det_u,det_v
		  alpha_u,alpha_v
		*/

		vector<double> tr_u(bases_num);
		vector<double> tr_v(bases_num);
		vector<double> det_u(bases_num);
		vector<double> det_v(bases_num);
		vector<double> alpha_u(bases_num);
		vector<double> alpha_v(bases_num);

		for (size_t j = 0; j < bases_num; j++)
		{
			tr_u[j] = 2 * (a*coeff_ac[j] + b * coeff_bd[j]);
			tr_v[j] = 2 * (c*coeff_ac[j] + d * coeff_bd[j]);
			det_u[j] = d * coeff_ac[j] - c * coeff_bd[j];
			det_v[j] = a * coeff_bd[j] - b * coeff_ac[j];
			alpha_u[j] = tr_u[j] + 2 * det_u[j];
			alpha_v[j] = tr_v[j] + 2 * det_v[j];
		}

		double tr = a * a + b * b + c * c + d * d;
		double det = a * d - b * c;
		double alpha = sqrt(tr + 2 * det);
		double beta = sqrt((a - d)*(a - d) + (b + c)*(b + c));
		double sig1 = (alpha + beta) / 2.0;
		double sig2 = (alpha - beta) / 2.0;
		double det_2 = 1.0 / det / det;
		double det_3 = det_2 / det;
		double det_4 = det_2 * det_2;

		//calc gradient
		for (size_t j = 0; j < bases_num; j++)
		{
			double g_u_ = f_area * (tr_u[j] * (1.0 + det_2) - 2.0*tr*det_u[j] * det_3);
			double g_v_ = f_area * (tr_v[j] * (1.0 + det_2) - 2.0*tr*det_v[j] * det_3);
			grad_compents_record.emplace_back(ids_bases[j], g_u_);
			grad_compents_record.emplace_back(ids_bases[j] + v_free_N, g_v_);
			//gradient_[ids_bases[j].first] -= g_u_;
			//gradient_[ids_bases[j].first + VAR_N] -= g_v_;
		}

		//calc hessian
		for (size_t ui = 0; ui < bases_num; ui++)
		{
			for (size_t uj = ui; uj < bases_num; uj++)
			{
				double tmp_1 = 2.0 * (coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj])*(1.0 + det_2);
				double h_ui_uj_plus = (tmp_1 - 2.0*det_3*(tr_u[ui] * det_u[uj] + tr_u[uj] * det_u[ui]) + 6.0*tr*det_4*det_u[ui] * det_u[uj]);
				double h_vi_vj_plus = (tmp_1 - 2.0*det_3*(tr_v[ui] * det_v[uj] + tr_v[uj] * det_v[ui]) + 6.0*tr*det_4*det_v[ui] * det_v[uj]);

				//ps_a[find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first)] += f_area * h_ui_uj_plus;
				//ps_a[find_idx_in_hess.coeff(ids_bases[ui].first + VAR_N, ids_bases[uj].first + VAR_N)] += f_area * h_vi_vj_plus;

				hess_compents_record.emplace_back(SM_find_id.coeff(ids_bases[ui], ids_bases[uj]), f_area*h_ui_uj_plus);
				hess_compents_record.emplace_back(SM_find_id.coeff(ids_bases[ui] + v_free_N, ids_bases[uj] + v_free_N), f_area*h_vi_vj_plus);

			}
			for (size_t vj = 0; vj < bases_num; vj++)
			{
				double h_ui_vj_plus = -2.0*(det_3*tr_u[ui] * det_v[vj] + tr * det_3*(coeff_bd[vj] * coeff_ac[ui] - coeff_ac[vj] * coeff_bd[ui]) + det_3 * tr_v[vj] * det_u[ui] - 3.0*tr*det_4*det_u[ui] * det_v[vj]);
				//ps_a[find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[vj].first + VAR_N)] += f_area * h_ui_vj_plus;
				hess_compents_record.emplace_back(SM_find_id.coeff(ids_bases[ui], ids_bases[vj] + v_free_N), f_area * h_ui_vj_plus);

			}
		}

		double sig_flag = 2.0*(sig1 + sig2 - 1.0 / sig1 / sig1 / sig1 - 1.0 / sig2 / sig2 / sig2);

		if (sig_flag < 0)
		{

			double alpha_1 = 1.0 / alpha;
			double alpha_3 = alpha_1 * alpha_1*alpha_1;
			for (size_t ui = 0; ui < bases_num; ui++)
			{
				for (size_t uj = ui; uj < bases_num; uj++)
				{
					double tmp_1 = (coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj])*alpha_1 / 2.0;
					double alpha_ui_uj_minus = sig_flag * (tmp_1 - alpha_u[ui] * alpha_u[uj] * alpha_3 / 8.0);
					double alpha_vi_vj_minus = sig_flag * (tmp_1 - alpha_v[ui] * alpha_v[uj] * alpha_3 / 8.0);
					//ps_a[find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[uj].first)] -= f_area * alpha_ui_uj_minus;
					//ps_a[find_idx_in_hess.coeff(ids_bases[ui].first + VAR_N, ids_bases[uj].first + VAR_N)] -= f_area * alpha_vi_vj_minus;
					hess_compents_record.emplace_back(SM_find_id.coeff(ids_bases[ui], ids_bases[uj]), -f_area * alpha_ui_uj_minus);
					hess_compents_record.emplace_back(SM_find_id.coeff(ids_bases[ui] + v_free_N, ids_bases[uj] + v_free_N), -f_area * alpha_vi_vj_minus);

				}
				for (size_t vj = 0; vj < bases_num; vj++)
				{
					double alpha_ui_vj_minus = sig_flag * ((coeff_bd[vj] * coeff_ac[ui] - coeff_ac[vj] * coeff_bd[ui])*alpha_1 / 2.0 - alpha_u[ui] * alpha_v[vj] * alpha_3 / 8.0);
					//ps_a[find_idx_in_hess.coeff(ids_bases[ui].first, ids_bases[vj].first + VAR_N)] -= f_area * alpha_ui_vj_minus;
					hess_compents_record.emplace_back(SM_find_id.coeff(ids_bases[ui], ids_bases[vj] + v_free_N), -f_area * alpha_ui_vj_minus);

				}
			}
		}

#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_local.a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				ps_local.rhs[var.first] -= var.second;
		}
	}

}

double LocalPara::max_step_general(const Eigen::MatrixXd & decrease_)
{
	int f0, f1, f2;
	double x02, y02, x12, y12, dx02, dy02, dx12, dy12;
	double a, b, c;
	double t_max = numeric_limits<double>::infinity();
	double t_tmp;
	for (int i = 0; i < F_p.rows(); i++)
	{
		f0 = F_p(i, 0);
		f1 = F_p(i, 1);
		f2 = F_p(i, 2);
		x02 = V_p(f0, 0) - V_p(f2, 0);
		y02 = V_p(f0, 1) - V_p(f2, 1);
		x12 = V_p(f1, 0) - V_p(f2, 0);
		y12 = V_p(f1, 1) - V_p(f2, 1);
		dx02 = decrease_(f0,0) - decrease_(f2,0);
		dy02 = decrease_(f0,1) - decrease_(f2,1);
		dx12 = decrease_(f1,0) - decrease_(f2,0);
		dy12 = decrease_(f1,1) - decrease_(f2,1);

		a = dx02 * dy12 - dx12 * dy02;
		b = x02 * dy12 + dx02 * y12 - y02 * dx12 - x12 * dy02;
		c = x02 * y12 - x12 * y02;
		t_tmp = get_smallest_pos_quad_zero(a, b, c);
		if (t_max > t_tmp)
		{
			t_max = t_tmp;
		}
	}
	return t_max;
}

void LocalPara::backtracking_general(const Eigen::MatrixXd & neg_grad_, const Eigen::MatrixXd & decrease_, double & step_)
{
	double h_ = 0.5;
	double c_ = 0.2;
	double tt = 0.0;

	tt = -(neg_grad_.col(0).dot(decrease_.col(0)));
	tt -= (neg_grad_.col(1).dot(decrease_.col(1)));
	//printf("tt==: %f\n", tt);

	double ex = calc_energy(V_p);
	//printf("ex==: %f\n", ex);

	Eigen::MatrixXd pos_new = V_p + step_ * decrease_;
	double e = calc_energy(pos_new);
	//printf("e==: %f\n", e);

	while (e > ex + step_ * c_*tt)
	{
		step_ *= h_;
		pos_new = V_p + step_ * decrease_;
		e = calc_energy(pos_new);
		//printf("e==: %f\n", e);
	}
}

double LocalPara::calc_energy(const Eigen::MatrixXd & pos_cur)
{
	double e_sum = 0.0;
	int f0, f1, f2;
	double det, tr, E_d;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	for (int i = 0; i < F_p.rows(); i++)
	{
		f0 = F_p(i, 0);
		f1 = F_p(i, 1);
		f2 = F_p(i, 2);
		q00 = pos_cur(f0,0) - pos_cur(f2,0);
		q01 = pos_cur(f1,0) - pos_cur(f2,0);
		q10 = pos_cur(f0,1) - pos_cur(f2,1);
		q11 = pos_cur(f1,1) - pos_cur(f2,1);
		p00 = src_inv_p[i][0]; p01 = src_inv_p[i][1]; p10 = src_inv_p[i][2]; p11 = src_inv_p[i][3];
		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;
		det = j00 * j11 - j01 * j10;
		tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
		E_d = tr + tr / (det*det);
		e_sum += weight_p[i] * E_d;
	}
	return e_sum;
}

void LocalPara::global2local_with_fixed_boundary()
{
	vector<vector<int>> boundaryloop;
	boundary_loop(F_src_p, boundaryloop);

	set<int> boundary_set;
	for (const auto&var : boundaryloop)
		boundary_set.insert(var.begin(), var.end());

	Eigen::VectorXi mark;
	mark.setZero(v_src_N);
	int shanchu = 0;
//	set<int> sanchu;
	for (size_t i = 0; i < F_src_p.rows(); i++)
	{
		bool is_free = false;
		for (size_t j = 0; j < F_src_p.cols(); j++)
		{
			if (boundary_set.find(F_src_p(i, j)) == boundary_set.end())
				mark(F_src_p(i, j)) = 1;
			else
				is_free = true;
		}
		if (is_free&&fids_p[i] < sD.mf_num)
		{
			f_visit[fids_p[i]]++;
//			sanchu.insert(i);
			shanchu++;
		}
	}
//	cout << "shanchu: " << shanchu << endl;

	v_free_N = mark.sum();
	v_g2l.resize(v_src_N);
	v_l2g.resize(v_free_N + boundary_set.size());
	v_w_N = v_l2g.size();
	int counts = 0;

	for (size_t i = 0; i < mark.size(); i++)
	{
		if (mark(i) == 1)
		{
			v_g2l(i) = counts;
			v_l2g(counts) = i;
			counts++;
		}
		else
		{
			v_g2l(i) = -1;
		}
	}

	for (const auto&vid : boundary_set)
	{
		v_g2l[vid] = counts;
		v_l2g(counts) = vid;
		counts++;
	}

	F_p.resize(F_src_p.rows(), F_src_p.cols());
	for (size_t i = 0; i < F_src_p.rows(); i++)
	{
		for (size_t j = 0; j < F_src_p.cols(); j++)
		{
			F_p(i, j) = v_g2l(F_src_p(i, j));
		}
	}

	V_p.resize(v_w_N, 2);
	for (size_t i = 0; i < v_w_N; i++)
	{
		V_p.row(i) = sD.w_uv.row(v_l2g[i]);
	}


	//Eigen::MatrixXi F_p_b;
	//F_p_b.resize(sanchu.size(), 3);
	//int i_c_ = 0;
	//for (auto itera = sanchu.begin(); itera != sanchu.end(); itera++)
	//{
	//	F_p_b.row(i_c_) = F_p.row(*itera);
	//	i_c_++;
	//}
	//writeObj(V_p, F_p, "part.obj");
	//writeObj(V_p, F_p_b, "part_b.obj");
	//system("pause");

}

void LocalPara::init()
{
	F_src_p.resize(fids_p.size(), 3);
	src_inv_p.resize(fids_p.size());
	weight_p.resize(fids_p.size());

	int fid_count = 0;
	weight_total = 0.0;
	for (auto fid_ = fids_p.begin(); fid_ != fids_p.end(); fid_++)
	{
		if (*fid_ < sD.mf_num)
		{
			F_src_p.row(fid_count) = sD.F.row(*fid_);
			weight_p[fid_count] = sD.src_weight[*fid_];
			src_inv_p[fid_count] = Eigen::Vector4d(sD.upd_inv00[*fid_], sD.upd_inv01[*fid_], sD.upd_inv10[*fid_], sD.upd_inv11[*fid_]);
		}
		else
		{
			F_src_p.row(fid_count) = sD.s_T.row(*fid_-sD.mf_num);
			weight_p[fid_count] = sD.scaf_weight;
			src_inv_p[fid_count] = sD.scaf_jac[*fid_ - sD.mf_num];
		}
		weight_total += weight_p[fid_count];
		fid_count++;
	}

	global2local_with_fixed_boundary();

	printf("local part == var_N: %d, beginning V_N: %d F_N: %d",v_free_N, v_w_N, fids_p.size());
}

void LocalPara::run()
{
	clock_t t1 = clock();
	init();
	if (v_free_N < 20)
		return;
	prepare_pardiso_with_fixed_boundary();
	ps_local.pardiso_init();
	double cur_conv_rate = 1.0;
	int cur_iter = 0;
	double e_cur = calc_energy(V_p);
	double e_init = e_cur;
	double e_pre;
	for (int iter = 0; iter < MAX_ITERATION_N; iter++)
	{
		e_pre = e_cur;
		solve_equation();
		e_cur = calc_energy(V_p);
		cur_iter++;
		cur_conv_rate = (e_pre - e_cur) / e_pre;
		if (cur_conv_rate < CONVERGENCE_RATE)
			break;
	}	
	update_src_pos();
	clock_t t2 = clock();
	printf(" == finish; time: %f , e_init: %f, e_end: %f \n", (t2 - t1) / 1000.0,e_init/weight_total,e_cur/weight_total);

}

void LocalPara::solve_equation()
{
	assembly_pardiso_with_fixed_boundary();
	ps_local.factorize();
	ps_local.pardiso_solver();
	Eigen::MatrixXd decrease_, neg_grad_;
	decrease_.setZero(v_w_N, 2);
	neg_grad_.setZero(v_w_N, 2);
	for (size_t i = 0; i < v_free_N; i++)
	{
		decrease_(i,0) = ps_local.result[i];
		decrease_(i,1) = ps_local.result[i + v_free_N];
		neg_grad_(i,0) = ps_local.rhs[i];
		neg_grad_(i,1) = ps_local.rhs[i + v_free_N];
	}
	double max_step_ = max_step_general(decrease_);
	double step_ = 0.95*max_step_;
	backtracking_general(neg_grad_, decrease_, step_);
	V_p += step_ * decrease_;
}

void LocalPara::update_src_pos()
{
	for (int i = 0; i < v_w_N; i++)
	{
		sD.w_uv.row(v_l2g[i]) = V_p.row(i);
		//V_src[v_l2g[i]] = V_p[i];
		//V_src[v_l2g[i] + v_src_N] = V_p[i + v_w_N];
	}

}

void LocalPara::prepare_pardiso_with_fixed_boundary()
{
	ps_local.ia.clear(); ps_local.ia.reserve(2 * v_free_N + 1);
	ps_local.ja.clear(); ps_local.ja.reserve(8 * v_free_N);

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(v_free_N);
	for (size_t i = 0; i < F_p.rows(); i++)
	{
		vector<int> vid;
		for (size_t j = 0; j < F_p.cols(); j++)
		{
			if (F_p(i, j) < v_free_N)
				vid.push_back(F_p(i, j));
		}
		for (size_t ii = 0; ii < vid.size(); ii++)
			for (size_t jj = 0; jj < vid.size(); jj++)
				VV_tmp[vid[ii]].insert(vid[jj]);
	}
	vector<Eigen::Triplet<int>> tripletlist;
	tripletlist.reserve(36 * v_free_N);

	for (int i = 0; i < v_free_N; i++)
	{
		ps_local.ia.push_back(ps_local.ja.size());
		int cur_id_start = ps_local.ia.back();
		vector<int> row_id;
		for (auto&var : VV_tmp[i])
		{
			row_id.push_back(var);
		}

		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		ps_local.ja.push_back(i);
		tripletlist.emplace_back(i, i, cur_id_start + dd);
		dd++;
		for (int k = std::distance(row_id.begin(), iter)+1; k < row_id.size(); k++)
		{
			ps_local.ja.push_back(row_id[k]);
			tripletlist.emplace_back(i, row_id[k], cur_id_start + dd);
			tripletlist.emplace_back(row_id[k], i, cur_id_start + dd);
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			ps_local.ja.push_back(row_id[k] + v_free_N);
			tripletlist.emplace_back(i, row_id[k] + v_free_N, cur_id_start + dd);
			tripletlist.emplace_back(row_id[k] + v_free_N, i, cur_id_start + dd);
			++dd;
		}
	}
	for (int i = v_free_N; i < 2 * v_free_N; i++)
	{
		ps_local.ia.push_back(ps_local.ja.size());
		int cur_id_start = ps_local.ia.back();
		vector<int> row_id;
		for (auto&var : VV_tmp[i - v_free_N])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i - v_free_N);

		int dd = 0;
		ps_local.ja.push_back(i);
		tripletlist.emplace_back(i, i, cur_id_start + dd);
		dd++;

		for (int k = std::distance(row_id.begin(), iter)+1; k < row_id.size(); k++)
		{
			ps_local.ja.push_back(row_id[k] + v_free_N);
			tripletlist.emplace_back(i, row_id[k] + v_free_N, cur_id_start + dd);
			tripletlist.emplace_back(row_id[k] + v_free_N, i, cur_id_start + dd);
			++dd;
		}
	}
	ps_local.ia.push_back(ps_local.ja.size());
	SM_find_id.resize(2 * v_free_N, 2 * v_free_N);
	SM_find_id.setFromTriplets(tripletlist.begin(), tripletlist.end());
	ps_local.nnz = ps_local.ja.size();
	ps_local.num = 2 * v_free_N;

}

