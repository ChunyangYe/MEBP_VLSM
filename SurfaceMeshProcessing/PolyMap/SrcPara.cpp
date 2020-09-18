#include "SrcPara.h"


SrcPara::SrcPara(SharedData & shardd_) :sD(shardd_)
{
}

SrcPara::~SrcPara()
{
	if (ps_local != NULL)
	{
		delete ps_local;
		ps_local = NULL;
	}
}

void SrcPara::assembly_pardiso_with_fixed_boundary()
{
	ps_local->a.clear();
	ps_local->a.assign(ps_local->ja.size(), 0.0);
	ps_local->rhs.clear();
	ps_local->rhs.assign(2 * v_free_N, 0.0);

#pragma omp parallel for
	for (int trii = 0; trii < sD.mf_num; trii++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		hess_compents_record.reserve(21);
		grad_compents_record.reserve(6);

		int f0 = sD.F(trii, 0);
		int f1 = sD.F(trii, 1);
		int f2 = sD.F(trii, 2);

		double f_area = sD.src_weight[trii];

		double p00 = sD.upd_inv00[trii];
		double p01 = sD.upd_inv01[trii];
		double p10 = sD.upd_inv10[trii];
		double p11 = sD.upd_inv11[trii];

		vector<double> coeff_ac;
		vector<double> coeff_bd;
		vector<int>ids_bases;
		coeff_ac.push_back(p00);
		coeff_bd.push_back(p01);
		ids_bases.push_back(f0);
		coeff_ac.push_back(p10);
		coeff_bd.push_back(p11);
		ids_bases.push_back(f1);
		coeff_ac.push_back(-(p00 + p10));
		coeff_bd.push_back(-(p01 + p11));
		ids_bases.push_back(f2);
		int bases_num = ids_bases.size();

		double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
		double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
		double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);
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
			for (size_t uj = 0; uj < bases_num; uj++)
			{
				if (ids_bases[uj] >= ids_bases[ui])
				{
					double tmp_1 = 2.0 * (coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj])*(1.0 + det_2);
					double h_ui_uj_plus = (tmp_1 - 2.0*det_3*(tr_u[ui] * det_u[uj] + tr_u[uj] * det_u[ui]) + 6.0*tr*det_4*det_u[ui] * det_u[uj]);
					double h_vi_vj_plus = (tmp_1 - 2.0*det_3*(tr_v[ui] * det_v[uj] + tr_v[uj] * det_v[ui]) + 6.0*tr*det_4*det_v[ui] * det_v[uj]);
					int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[uj]);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, f_area*h_ui_uj_plus);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui] + v_free_N] + dd_, f_area*h_vi_vj_plus);
				}
			}
			for (size_t vj = 0; vj < bases_num; vj++)
			{
				double h_ui_vj_plus = -2.0*(det_3*tr_u[ui] * det_v[vj] + tr * det_3*(coeff_bd[vj] * coeff_ac[ui] - coeff_ac[vj] * coeff_bd[ui]) + det_3 * tr_v[vj] * det_u[ui] - 3.0*tr*det_4*det_u[ui] * det_v[vj]);

				int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[vj] + v_free_N);
				hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, f_area * h_ui_vj_plus);

			}
		}

		double sig_flag = 2.0*(sig1 + sig2 - 1.0 / sig1 / sig1 / sig1 - 1.0 / sig2 / sig2 / sig2);

		if (sig_flag < 0)
		{

			double alpha_1 = 1.0 / alpha;
			double alpha_3 = alpha_1 * alpha_1*alpha_1;
			for (size_t ui = 0; ui < bases_num; ui++)
			{
				for (size_t uj = 0; uj < bases_num; uj++)
				{
					if (ids_bases[uj] >= ids_bases[ui])
					{
						double tmp_1 = (coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj])*alpha_1 / 2.0;
						double alpha_ui_uj_minus = sig_flag * (tmp_1 - alpha_u[ui] * alpha_u[uj] * alpha_3 / 8.0);
						double alpha_vi_vj_minus = sig_flag * (tmp_1 - alpha_v[ui] * alpha_v[uj] * alpha_3 / 8.0);

						int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[uj]);
						hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, -f_area * alpha_ui_uj_minus);
						hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui] + v_free_N] + dd_, -f_area * alpha_vi_vj_minus);
					}

				}
				for (size_t vj = 0; vj < bases_num; vj++)
				{
					double alpha_ui_vj_minus = sig_flag * ((coeff_bd[vj] * coeff_ac[ui] - coeff_ac[vj] * coeff_bd[ui])*alpha_1 / 2.0 - alpha_u[ui] * alpha_v[vj] * alpha_3 / 8.0);

					int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[vj] + v_free_N);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, -f_area * alpha_ui_vj_minus);
				}
			}
		}

#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_local->a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				ps_local->rhs[var.first] -= var.second;
		}
	}
	assembly_pardiso_scaffold();
}

void SrcPara::assembly_pardiso_scaffold()
{
#pragma omp parallel for
	for (int trii = 0; trii < sD.sf_num; trii++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		hess_compents_record.reserve(21);
		grad_compents_record.reserve(6);

		int f0 = sD.s_T(trii, 0);
		int f1 = sD.s_T(trii, 1);
		int f2 = sD.s_T(trii, 2);

		double f_area = sD.scaf_weight;

		double p00 = sD.scaf_jac[trii][0];
		double p01 = sD.scaf_jac[trii][1];
		double p10 = sD.scaf_jac[trii][2];
		double p11 = sD.scaf_jac[trii][3];

		double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
		double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
		double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);

		vector<double> coeff_ac;
		vector<double> coeff_bd;
		vector<int>ids_bases;
		if (src2var[f0] != -1)
		{
			coeff_ac.push_back(p00);
			coeff_bd.push_back(p01);
			ids_bases.push_back(src2var[f0]);
		}
		if (src2var[f1] != -1)
		{
			coeff_ac.push_back(p10);
			coeff_bd.push_back(p11);
			ids_bases.push_back(src2var[f1]);
		}
		if (src2var[f2] != -1)
		{
			coeff_ac.push_back(-(p00 + p10));
			coeff_bd.push_back(-(p01 + p11));
			ids_bases.push_back(src2var[f2]);
		}
		int bases_num = ids_bases.size();

		if (bases_num == 0)
			continue;


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
			for (size_t uj = 0; uj < bases_num; uj++)
			{
				if (ids_bases[uj] >= ids_bases[ui])
				{
					double tmp_1 = 2.0 * (coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj])*(1.0 + det_2);
					double h_ui_uj_plus = (tmp_1 - 2.0*det_3*(tr_u[ui] * det_u[uj] + tr_u[uj] * det_u[ui]) + 6.0*tr*det_4*det_u[ui] * det_u[uj]);
					double h_vi_vj_plus = (tmp_1 - 2.0*det_3*(tr_v[ui] * det_v[uj] + tr_v[uj] * det_v[ui]) + 6.0*tr*det_4*det_v[ui] * det_v[uj]);


					int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[uj]);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, f_area*h_ui_uj_plus);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui] + v_free_N] + dd_, f_area*h_vi_vj_plus);
				}
			}
			for (size_t vj = 0; vj < bases_num; vj++)
			{
				double h_ui_vj_plus = -2.0*(det_3*tr_u[ui] * det_v[vj] + tr * det_3*(coeff_bd[vj] * coeff_ac[ui] - coeff_ac[vj] * coeff_bd[ui]) + det_3 * tr_v[vj] * det_u[ui] - 3.0*tr*det_4*det_u[ui] * det_v[vj]);

				int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[vj] + v_free_N);
				hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, f_area * h_ui_vj_plus);

			}
		}

		double sig_flag = 2.0*(sig1 + sig2 - 1.0 / sig1 / sig1 / sig1 - 1.0 / sig2 / sig2 / sig2);

		if (sig_flag < 0)
		{
			double alpha_1 = 1.0 / alpha;
			double alpha_3 = alpha_1 * alpha_1*alpha_1;
			for (size_t ui = 0; ui < bases_num; ui++)
			{
				for (size_t uj = 0; uj < bases_num; uj++)
				{
					if (ids_bases[uj] >= ids_bases[ui])
					{
						double tmp_1 = (coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj])*alpha_1 / 2.0;
						double alpha_ui_uj_minus = sig_flag * (tmp_1 - alpha_u[ui] * alpha_u[uj] * alpha_3 / 8.0);
						double alpha_vi_vj_minus = sig_flag * (tmp_1 - alpha_v[ui] * alpha_v[uj] * alpha_3 / 8.0);

						int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[uj]);
						hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, -f_area * alpha_ui_uj_minus);
						hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui] + v_free_N] + dd_, -f_area * alpha_vi_vj_minus);
					}

				}
				for (size_t vj = 0; vj < bases_num; vj++)
				{
					double alpha_ui_vj_minus = sig_flag * ((coeff_bd[vj] * coeff_ac[ui] - coeff_ac[vj] * coeff_bd[ui])*alpha_1 / 2.0 - alpha_u[ui] * alpha_v[vj] * alpha_3 / 8.0);

					int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[vj] + v_free_N);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, -f_area * alpha_ui_vj_minus);
				}
			}
		}


#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_local->a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				ps_local->rhs[var.first] -= var.second;
		}
	}

}

void SrcPara::assemble_hessian_gradient_slim()
{
	ps_local->a.clear();
	ps_local->a.assign(ps_local->ja.size(), 0.0);
	ps_local->rhs.clear();
	ps_local->rhs.assign(2 * v_free_N, 0.0);
	int total_f_num = sD.mf_num + sD.s_T.rows();
#pragma omp parallel for
	for (int trii = 0; trii < total_f_num; trii++)
	{
		vector<pair<int, double>> hess_compents_record;
		vector<pair<int, double>> grad_compents_record;

		hess_compents_record.reserve(21);
		grad_compents_record.reserve(6);

		int f0, f1, f2;
		double p00, p01, p10, p11, f_area;
		vector<double> coeff_ac;
		vector<double> coeff_bd;
		vector<int>ids_bases;

		if (trii < sD.mf_num)
		{
			f0 = sD.F(trii, 0);
			f1 = sD.F(trii, 1);
			f2 = sD.F(trii, 2);
			f_area = sD.src_weight[trii];

			p00 = sD.upd_inv00[trii];
			p01 = sD.upd_inv01[trii];
			p10 = sD.upd_inv10[trii];
			p11 = sD.upd_inv11[trii];

			coeff_ac.push_back(p00);
			coeff_bd.push_back(p01);
			ids_bases.push_back(f0);
			coeff_ac.push_back(p10);
			coeff_bd.push_back(p11);
			ids_bases.push_back(f1);
			coeff_ac.push_back(-(p00 + p10));
			coeff_bd.push_back(-(p01 + p11));
			ids_bases.push_back(f2);
		}
		else
		{
			int trij = trii - sD.mf_num;
			f0 = sD.s_T(trij, 0);
			f1 = sD.s_T(trij, 1);
			f2 = sD.s_T(trij, 2);
			f_area = sD.scaf_weight;
			p00 = sD.scaf_jac[trij][0];
			p01 = sD.scaf_jac[trij][1];
			p10 = sD.scaf_jac[trij][2];
			p11 = sD.scaf_jac[trij][3];

			if (src2var[f0] != -1)
			{
				coeff_ac.push_back(p00);
				coeff_bd.push_back(p01);
				ids_bases.push_back(src2var[f0]);
			}
			if (src2var[f1] != -1)
			{
				coeff_ac.push_back(p10);
				coeff_bd.push_back(p11);
				ids_bases.push_back(src2var[f1]);
			}
			if (src2var[f2] != -1)
			{
				coeff_ac.push_back(-(p00 + p10));
				coeff_bd.push_back(-(p01 + p11));
				ids_bases.push_back(src2var[f2]);
			}
		}

		int bases_num = ids_bases.size();
		if (bases_num == 0)
			continue;

		double q00 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
		double q10 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
		double q01 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		double q11 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);

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

		double tr = a * a + b * b + c * c + d * d;
		double det = a * d - b * c;
		double alpha = sqrt(tr + 2 * det);
		double beta = sqrt((a - d)*(a - d) + (b + c)*(b + c));
		double sig1 = (alpha + beta) / 2.0;
		double sig2 = (alpha - beta) / 2.0;
		double det_1 = 1.0 / det;
		double det_2plus1 = det_1 * det_1;
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
			double g_u_ = r0 * coeff_ac[j] + r1 * coeff_bd[j];
			double g_v_ = r2 * coeff_ac[j] + r3 * coeff_bd[j];
			grad_compents_record.emplace_back(ids_bases[j], g_u_);
			grad_compents_record.emplace_back(ids_bases[j] + v_free_N, g_v_);
		}

		//calc hessian
		for (size_t ui = 0; ui < bases_num; ui++)
		{
			for (size_t uj = 0; uj < bases_num; uj++)
			{
				double h_ui_uj_plus = coeff_ac[ui] * coeff_ac[uj] + coeff_bd[ui] * coeff_bd[uj];

				hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + SM_find_id.coeff(ids_bases[ui], ids_bases[uj] + v_free_N), w2*h_ui_uj_plus);

				if (ids_bases[uj] >= ids_bases[ui])
				{
					int dd_ = SM_find_id.coeff(ids_bases[ui], ids_bases[uj]);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui]] + dd_, w1*h_ui_uj_plus);
					hess_compents_record.emplace_back(ps_local->ia[ids_bases[ui] + v_free_N] + dd_, w3*h_ui_uj_plus);
				}
			}
		}

#pragma omp critical
		{
			for (const auto&var : hess_compents_record)
				ps_local->a[var.first] += var.second;
			for (const auto&var : grad_compents_record)
				ps_local->rhs[var.first] -= var.second;
		}
	}

}

double SrcPara::max_step_general(const Eigen::MatrixXd & decrease_)
{
	int f0, f1, f2;
	double x02, y02, x12, y12, dx02, dy02, dx12, dy12;
	double a, b, c;
	double t_max = numeric_limits<double>::infinity();
	double t_tmp;
	for (int i = 0; i < sD.mf_num; i++)
	{
		f0 = sD.F(i, 0);
		f1 = sD.F(i, 1);
		f2 = sD.F(i, 2);
		x02 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
		y02 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
		x12 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		y12 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);
		dx02 = decrease_(f0, 0) - decrease_(f2, 0);
		dy02 = decrease_(f0, 1) - decrease_(f2, 1);
		dx12 = decrease_(f1, 0) - decrease_(f2, 0);
		dy12 = decrease_(f1, 1) - decrease_(f2, 1);

		a = dx02 * dy12 - dx12 * dy02;
		b = x02 * dy12 + dx02 * y12 - y02 * dx12 - x12 * dy02;
		c = x02 * y12 - x12 * y02;
		t_tmp = get_smallest_pos_quad_zero(a, b, c);
		if (t_max > t_tmp)
		{
			t_max = t_tmp;
		}
	}
	for (int i = 0; i < sD.sf_num; i++)
	{
		f0 = sD.s_T(i, 0);
		f1 = sD.s_T(i, 1);
		f2 = sD.s_T(i, 2);
		x02 = sD.w_uv(f0, 0) - sD.w_uv(f2, 0);
		y02 = sD.w_uv(f0, 1) - sD.w_uv(f2, 1);
		x12 = sD.w_uv(f1, 0) - sD.w_uv(f2, 0);
		y12 = sD.w_uv(f1, 1) - sD.w_uv(f2, 1);
		dx02 = decrease_(f0, 0) - decrease_(f2, 0);
		dy02 = decrease_(f0, 1) - decrease_(f2, 1);
		dx12 = decrease_(f1, 0) - decrease_(f2, 0);
		dy12 = decrease_(f1, 1) - decrease_(f2, 1);

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

void SrcPara::backtracking_general(const Eigen::MatrixXd & neg_grad_, const Eigen::MatrixXd & decrease_, double & step_)
{
	double h_ = 0.5;
	double c_ = 0.2;
	double tt = 0.0;

	tt = -(neg_grad_.col(0).dot(decrease_.col(0)));
	tt -= (neg_grad_.col(1).dot(decrease_.col(1)));

	if (tt > 0)
	{
		step_ = 0.0;
		return;
	}

	double ex = calc_energy(sD.w_uv);
	//printf("ex==: %f\n", ex);

	Eigen::MatrixXd pos_new = sD.w_uv + step_ * decrease_;
	double e = calc_energy(pos_new);
	//printf("e==: %f\n", e);

	while (e > ex + step_ * c_*tt)
	{
		step_ *= h_;
		pos_new = sD.w_uv + step_ * decrease_;
		e = calc_energy(pos_new);
		//printf("e==: %f\n", e);
	}
}

double SrcPara::calc_energy(const Eigen::MatrixXd & pos_cur, bool is_add_scaf)
{
	double e_sum = 0.0;

#pragma omp parallel for
	for (int i = 0; i < sD.mf_num; i++)
	{
		int f0 = sD.F(i, 0);
		int f1 = sD.F(i, 1);
		int f2 = sD.F(i, 2);
		double q00 = pos_cur(f0, 0) - pos_cur(f2, 0);
		double q01 = pos_cur(f1, 0) - pos_cur(f2, 0);
		double q10 = pos_cur(f0, 1) - pos_cur(f2, 1);
		double q11 = pos_cur(f1, 1) - pos_cur(f2, 1);
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
		double E_d = tr + tr / (det*det);
#pragma omp critical
		{
			e_sum += sD.src_weight[i] * E_d;
		}
	}

	if (is_add_scaf)
	{
#pragma omp parallel for
		for (int i = 0; i < sD.sf_num; i++)
		{
			int f0 = sD.s_T(i, 0);
			int f1 = sD.s_T(i, 1);
			int f2 = sD.s_T(i, 2);
			double q00 = pos_cur(f0, 0) - pos_cur(f2, 0);
			double q01 = pos_cur(f1, 0) - pos_cur(f2, 0);
			double q10 = pos_cur(f0, 1) - pos_cur(f2, 1);
			double q11 = pos_cur(f1, 1) - pos_cur(f2, 1);
			double p00 = sD.scaf_jac[i][0];
			double p01 = sD.scaf_jac[i][1];
			double p10 = sD.scaf_jac[i][2];
			double p11 = sD.scaf_jac[i][3];
			double j00 = p00 * q00 + p10 * q01;
			double j01 = p01 * q00 + p11 * q01;
			double j10 = p00 * q10 + p10 * q11;
			double j11 = p01 * q10 + p11 * q11;
			double det = j00 * j11 - j01 * j10;
			double tr = j00 * j00 + j01 * j01 + j10 * j10 + j11 * j11;
			double E_d = tr + tr / (det*det);
#pragma omp critical
			{
				e_sum += sD.scaf_weight* E_d;
			}
		}
	}

	return e_sum;
}

void SrcPara::global2local_with_fixed_boundary()
{
	v_w_N = sD.V_num;
	int bnd_n = sD.frame_ids.size();
	v_free_N = v_w_N - bnd_n;
	var2src.resize(v_w_N - bnd_n);

#pragma omp parallel for
	for (int i = 0; i < sD.mv_num; i++)
	{
		var2src[i] = i;
	}
#pragma omp parallel for
	for (int j = sD.mv_num + bnd_n; j < v_w_N; j++)
	{
		var2src[j - bnd_n] = j;
	}

	src2var.resize(v_w_N, -1);
#pragma omp parallel for
	for (int i = 0; i < var2src.size(); i++)
	{
		src2var[var2src[i]] = i;
	}
	printf("v_w_N %d; bnd_n %d; v_free_N %d;\n", v_w_N, bnd_n, v_free_N);
}

void SrcPara::init()
{
//	ps_local = NULL;

	printf("run_src_begin==========\n");
	LOG(INFO) << "run_src_begin==========";
}

void SrcPara::run()
{
	if (ps_local != NULL)
	{
		delete ps_local;
		ps_local = NULL;
	}
	ps_local = new PardisoSolver();
	if (!is_pardiso_solver)
		ps_local->set_as_the_container();
	clock_t t1s = clock();
	init();
	double cur_conv_rate = 1.0;
	int cur_iter = 0;
	double e_cur = sD.Energy_cur;
	double e_pre;
	for (int iter = 0; iter < MAX_ITERATION_N; iter++)
	{
		clock_t t1 = clock();

		e_pre = e_cur;
		sD.mesh_improve(true);
		sD.adjust_scaf_weight(sD.Energy_cur*sD.mesh_measure / sD.sf_num / 100.0);

		sD.update_scaf_jacobian();
		global2local_with_fixed_boundary();

		solve_equation();
		sD.Energy_cur = calc_energy(sD.w_uv, false);
		e_cur = sD.Energy_cur;
		sD.energy_record.push_back(sD.Energy_cur);

		clock_t t2 = clock();

		cur_iter++;
		printf("iter %d ,== time: %f , energy: %f, e_decrease: %f\n", cur_iter, (t2 - t1) / 1000.0, e_cur, e_pre - e_cur);
		LOG(INFO) << "iter " << cur_iter << ",== time: " << (t2 - t1) / 1000.0 << " , energy: " << e_cur << ", e_decrease: " << e_pre - e_cur;
		sD.energy_record3.emplace_back(sD.Energy_cur, (t2 - t1) / 1000.0, 0);

		cur_conv_rate = (e_pre - e_cur) / e_pre;
		if (cur_conv_rate < CONVERGENCE_RATE)
			break;
	}
	clock_t t2s = clock();
	printf("run_src_finish; sum_time: %f ,energy_conv: %f\n", (t2s - t1s) / 1000.0, sD.Energy_cur);
	LOG(INFO) << "run_src_finish; sum_time: " << (t2s - t1s) / 1000.0 << " ,energy_conv: " << sD.Energy_cur;
}

void SrcPara::set_method(int meth_type)
{
	//out-of_core
	if (meth_type == 0)
	{
		printf("out-of_core begin-----------\n");
		LOG(INFO) << "out-of_core begin-----------";
		is_amgcl_solver = false;
		is_pardiso_solver = true;
	}
	else if(meth_type==1)//amgcl
	{
		printf("amgcl begin-----------\n");
		LOG(INFO) << "amgcl begin-----------";
		is_amgcl_solver = true;
		is_pardiso_solver = false;
	}
}

void SrcPara::solve_equation()
{
	prepare_pardiso_with_fixed_boundary();
	assembly_pardiso_with_fixed_boundary();
//	assemble_hessian_gradient_slim();
	if (is_pardiso_solver)
	{
		ps_local->pardiso_init();
		ps_local->factorize();
		ps_local->pardiso_solver();
	}
	else
	{
		if (is_amgcl_solver)
		{
			amgcl_solver();
		}
		else
		{
			GS_solver();
		}
		ps_local->result = result_gs;
		result_gs.clear();
	}

	Eigen::MatrixXd decrease_, neg_grad_;
	decrease_.setZero(v_w_N, 2);
	neg_grad_.setZero(v_w_N, 2);

#pragma omp parallel for
	for (int i = 0; i < v_free_N; i++)
	{
		decrease_(var2src[i], 0) = ps_local->result[i];
		decrease_(var2src[i], 1) = ps_local->result[i + v_free_N];
		neg_grad_(var2src[i], 0) = ps_local->rhs[i];
		neg_grad_(var2src[i], 1) = ps_local->rhs[i + v_free_N];
	}

	//ps_local->free_numerical_factorization_memory();
	if (is_pardiso_solver)
	{
		delete ps_local;
		ps_local = NULL;
		ps_local = new PardisoSolver();
	}


	double max_step_ = max_step_general(decrease_);
	double step_ = min(1.0, 0.95*max_step_);
	backtracking_general(neg_grad_, decrease_, step_);
	sD.w_uv += step_ * decrease_;
}

void SrcPara::graph_color_GS(const std::vector<std::set<int>>& VVs)
{
	printf("Color====================\n");
	phase.clear();
	int nsizes = v_free_N;
	std::vector<int> color(2 * nsizes, -1);
	std::vector<uint8_t>possible_color;

	uint8_t ncolors = 0u;
	for (int i = 0; i < nsizes; i++)
	{
		std::fill(possible_color.begin(), possible_color.end(), 1u);
		int c;
		for (auto itvv : VVs[i])
		{
			c = color[itvv];
			if (c >= 0)
			{
				possible_color[c] = 0;
			}
			c = color[itvv + nsizes];
			if (c >= 0)
			{
				possible_color[c] = 0;
			}
		}
		int color_id = -1;
		for (auto j = 0; j < possible_color.size(); j++)
		{
			if (possible_color[j] != 0)
			{
				color_id = j;
				break;
			}
		}
		if (color_id < 0)
		{
			color_id = ncolors++;
			possible_color.resize(ncolors);
		}
		color[i] = color_id;
	}
	for (int i = nsizes; i < 2 * nsizes; i++)
	{
		std::fill(possible_color.begin(), possible_color.end(), 1u);
		int c;
		for (auto itvv : VVs[i - nsizes])
		{
			c = color[itvv];
			if (c >= 0)
			{
				possible_color[c] = 0;
			}
			c = color[itvv + nsizes];
			if (c >= 0)
			{
				possible_color[c] = 0;
			}
		}
		int color_id = -1;
		for (auto j = 0; j < possible_color.size(); j++)
		{
			if (possible_color[j] != 0)
			{
				color_id = j;
				break;
			}
		}
		if (color_id < 0)
		{
			color_id = ncolors++;
			possible_color.resize(ncolors);
		}
		color[i] = color_id;
	}

	phase.resize(ncolors);
	for (int i = 0; i < 2 * nsizes; i++)
	{
		phase[color[i]].push_back(i);
	}

	printf("max_color: %d vs colors_size: %d\n", ncolors, phase.size());
	for (size_t i = 0; i < ncolors; i++)
	{
		cout << phase[i].size() << endl;

	}

}

void SrcPara::GS_solver()
{
	vector<vector<pair<int, int>>> mat_L(2 * v_free_N);
	for (int i = 0; i < 2 * v_free_N; i++)
	{
		for (int j = ps_local->ia[i] + 1; j < ps_local->ia[i + 1]; j++)
		{
			mat_L[ps_local->ja[j]].emplace_back(i, j);
		}
	}

	result_gs.clear();
	result_gs.resize(2 * v_free_N, 0.0);
	vector<double> result_prev;
	for (int iter_gs = 0; iter_gs < 5000; iter_gs++)
	{
		result_prev = result_gs;
		for (int c_idx = 0; c_idx < phase.size(); c_idx++)
		{
#pragma omp parallel for
			for (int i = 0; i < phase[c_idx].size(); i++)
			{
				int cur_id = phase[c_idx][i];
				double resi_ = ps_local->rhs[cur_id];
				for (int li = 0; li < mat_L[cur_id].size(); li++)
				{
					resi_ -= result_gs[mat_L[cur_id][li].first] * ps_local->a[mat_L[cur_id][li].second];
				}
				for (int ri = ps_local->ia[cur_id] + 1; ri < ps_local->ia[cur_id + 1]; ri++)
				{
					resi_ -= result_gs[ps_local->ja[ri]] * ps_local->a[ri];
				}
				result_gs[cur_id] = resi_ / ps_local->a[ps_local->ia[cur_id]];
			}
		}

		double resid_ = 0.0;
		double norm_ = 0.0;
		for (int i = 0; i < 2 * v_free_N; i++)
		{
			double tmp_ = result_gs[i] - result_prev[i];
			resid_ += tmp_ * tmp_;
			norm_ += result_gs[i] * result_gs[i];
		}
		double conv_ = sqrt(resid_ / norm_);
		if (iter_gs % 100 == 0)
			printf("iter_gs: %d, resid_: %f, norm_: %f, conv_rate: %f\n", iter_gs, resid_, norm_, conv_);
		if (conv_ < 1e-3)
		{
			printf("iter_gs: %d, resid_: %f, norm_: %f, conv_rate: %f=========\n", iter_gs, resid_, norm_, conv_);
			break;
		}
	}
	phase.clear();
	mat_L.clear();
	result_prev.clear();

}


void SrcPara::amg_write()
{
	//string fname = "mat";
	//std::ofstream f(fname.c_str());

	//// Banner
	//f << "%%MatrixMarket matrix coordinate real symmetric\n";
	//int rows = ps_local->num;
	//int cols = ps_local->num;
	//int nnz = ps_local->nnz;
	//// Sizes
	//f << rows << " " << cols << " " << nnz << "\n";

	//// Data
	//for (size_t i = 0; i < rows; ++i) {
	//	for (int j = ps_local->ia[i]; j < ps_local->ia[i + 1];j++) {			
	//		f << i + 1 << " " << ps_local->ja[j] + 1 << " " << ps_local->a[j] << "\n";
	//	}
	//}
	//f.close();

	string fname = "rhs";
	std::ofstream f(fname.c_str());

	f << "%%MatrixMarket matrix array real general\n";

	int rows = ps_local->num;

	// Sizes
	f << rows << " " << 1 << "\n";

	// Data
	for (size_t i = 0; i < rows; ++i) {
		f << ps_local->rhs[i] << "\n";
	}
	f.close();

}

void SrcPara::prepare_pardiso_with_fixed_boundary()
{

	ps_local->ia.clear(); ps_local->ia.reserve(2 * v_free_N + 1);
	ps_local->ja.clear(); ps_local->ja.reserve(8 * v_free_N);

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(v_free_N);
	for (int i = 0; i < sD.mf_num; i++)
	{
		vector<int> vid;
		for (size_t j = 0; j < 3; j++)
		{
			vid.push_back(sD.F(i, j));
		}
		for (size_t ii = 0; ii < vid.size(); ii++)
			for (size_t jj = 0; jj < vid.size(); jj++)
				VV_tmp[vid[ii]].insert(vid[jj]);
	}

	for (int i = 0; i < sD.sf_num; i++)
	{
		vector<int> vid;
		for (size_t j = 0; j < 3; j++)
		{
			int s_id = src2var[sD.s_T(i, j)];
			if (s_id != -1)
				vid.push_back(s_id);
		}
		for (size_t ii = 0; ii < vid.size(); ii++)
			for (size_t jj = 0; jj < vid.size(); jj++)
				VV_tmp[vid[ii]].insert(vid[jj]);
	}
	vector<Eigen::Triplet<int>> tripletlist;
	tripletlist.reserve(36 * v_free_N);

	for (int i = 0; i < v_free_N; i++)
	{
		ps_local->ia.push_back(ps_local->ja.size());
		vector<int> row_id(VV_tmp[i].begin(), VV_tmp[i].end());
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			ps_local->ja.push_back(row_id[k]);
			tripletlist.emplace_back(i, row_id[k], dd);
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			ps_local->ja.push_back(row_id[k] + v_free_N);
			tripletlist.emplace_back(i, row_id[k] + v_free_N, dd);
			++dd;
		}
	}
	for (int i = v_free_N; i < 2 * v_free_N; i++)
	{
		ps_local->ia.push_back(ps_local->ja.size());
		vector<int> row_id(VV_tmp[i - v_free_N].begin(), VV_tmp[i - v_free_N].end());
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i - v_free_N);

		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			ps_local->ja.push_back(row_id[k] + v_free_N);
		}
	}
	ps_local->ia.push_back(ps_local->ja.size());
	SM_find_id.resize(v_free_N, 2 * v_free_N);
	SM_find_id.setFromTriplets(tripletlist.begin(), tripletlist.end());
	ps_local->nnz = ps_local->ja.size();
	ps_local->num = 2 * v_free_N;
	if (!is_pardiso_solver && !is_amgcl_solver)
	{
		graph_color_GS(VV_tmp);
	}
	VV_tmp.clear();
	tripletlist.clear();

}

void SrcPara::amgcl_solver()
{
	int varn_ = 2 * v_free_N;
	int nnz_ = ps_local->a.size();
	vector<vector<pair<int, int>>> mat_L(2 * v_free_N);
	for (int i = 0; i < 2 * v_free_N; i++)
	{
		for (int j = ps_local->ia[i] + 1; j < ps_local->ia[i + 1]; j++)
		{
			mat_L[ps_local->ja[j]].emplace_back(i, j);
		}
	}

	result_gs.clear();
	result_gs.resize(2 * v_free_N, 0.0);

	vector<int> amg_ia(2 * v_free_N + 1, 0), amg_ja(nnz_ * 2 - varn_, 0);
	vector<double> amg_val(nnz_ * 2 - varn_);

	int count_ = 0;
	for (int cur_id = 0; cur_id < varn_; cur_id++)
	{
		int cur_id_nnz = mat_L[cur_id].size() + ps_local->ia[cur_id + 1] - ps_local->ia[cur_id];
		amg_ia[cur_id + 1] = amg_ia[cur_id] + cur_id_nnz;
		for (int li = 0; li < mat_L[cur_id].size(); li++)
		{
			amg_ja[count_] = mat_L[cur_id][li].first;
			amg_val[count_] = ps_local->a[mat_L[cur_id][li].second];
			count_++;
		}
		for (int ri = ps_local->ia[cur_id]; ri < ps_local->ia[cur_id + 1]; ri++)
		{
			amg_ja[count_] = ps_local->ja[ri];
			amg_val[count_] = ps_local->a[ri];
			count_++;
		}
	}
	//amg_write();
	//system("pause");
	mat_L.clear();
	ps_local->ia.clear();
	ps_local->ja.clear();
	ps_local->a.clear();
	using namespace amgcl;
	typedef amgcl::backend::builtin<double> Backend;
	boost::property_tree::ptree prm;
	amgcl::put(prm, "solver.tol=1e-4");
	amgcl::put(prm, "solver.type=bicgstab");
	//amgcl::put(prm, "solver.maxiter = 100");
	prm.put("solver.maxiter", 100);
	Backend::params bprm;
	typedef amgcl::make_solver<
		amgcl::runtime::preconditioner<Backend>,
		amgcl::runtime::solver::wrapper<Backend>
	> Solver;



	profiler<> prof;
	prof.tic("setup");
	Solver solve(std::tie(varn_, amg_ia, amg_ja, amg_val), prm, bprm);
	double tm_setup = prof.toc("setup");

	std::cout << solve << std::endl;

	LOG(INFO) << solve << std::endl;


	//auto f_b = Backend::copy_vector(ps_local->rhs, bprm);
	//auto x_b = Backend::copy_vector(result_gs, bprm);

	int iters;
	double error;
	prof.tic("solve");
	//	std::tie(iters, error) = solve(*f_b, *x_b);
	std::tie(iters, error) = solve(ps_local->rhs, result_gs);
	prof.toc("solve");

	std::cout
		<< "iters: " << iters << std::endl
		<< "error: " << error << std::endl
		<< prof << std::endl;

	LOG(INFO)
		<< "iters: " << iters << " error: " << error << std::endl
		<< prof << std::endl;

}
