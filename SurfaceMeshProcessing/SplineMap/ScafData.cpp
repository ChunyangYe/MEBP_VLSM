#include "ScafData.h"
#include <iostream>


using namespace std;
using namespace Eigen;

void ScafData::update_scaffold()
{
  mv_num = m_V.rows();
  mf_num = m_T.rows();

  v_num = w_uv.rows();
  sf_num = s_T.rows();

  sv_num = v_num - mv_num;
  f_num = sf_num + mf_num;

  s_M = Eigen::VectorXd::Constant(sf_num, scaffold_factor);
}

void ScafData::mesh_improve(bool in_packing = false)
{
	MatrixXd m_uv = w_uv.topRows(mv_num);
	MatrixXd V_bnd;
	V_bnd.resize(internal_bnd.size(), 2);
	for (int i = 0; i < internal_bnd.size(); i++) // redoing step 1.
	{
		V_bnd.row(i) = m_uv.row(internal_bnd(i));
	}

	if (rect_frame_V.size() == 0) {
		Matrix2d ob;// = rect_corners;
		VectorXd uv_max = m_uv.colwise().maxCoeff();
		VectorXd uv_min = m_uv.colwise().minCoeff();
		VectorXd uv_mid = (uv_max + uv_min) / 2.;
		double maxlen = (uv_max - uv_mid).maxCoeff();
		uv_min(0) = uv_mid(0) - maxlen;
		uv_min(1) = uv_mid(1) - maxlen;
		uv_max(0) = uv_mid(0) + maxlen;
		uv_max(1) = uv_mid(1) + maxlen;

		//double scaf_range = 3;
		Eigen::Array2d scaf_range(3, 3);
		ob.row(0) = uv_mid.array() + scaf_range * ((uv_min - uv_mid).array());
		ob.row(1) = uv_mid.array() + scaf_range * ((uv_max - uv_mid).array());

		ob.row(1) = ob.row(1) - ob.row(0);
		for (int i = 0; i < mv_num; ++i)
		{
			m_uv.row(i) = m_uv.row(i) - ob.row(0);
			w_uv.row(i) = w_uv.row(i) - ob.row(0);
		}
		ob.row(0).setZero();

		Vector2d rect_len;
		rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);
		int frame_points = 100;
		interval = rect_len(0);

		rect_frame_V.resize(4 * frame_points, 2);
		for (int i = 0; i < frame_points; i++)
		{
			// 0,0;0,1
			rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
			// 0,0;1,1
			rect_frame_V.row(i + frame_points)
				<< ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
			// 1,0;1,1
			rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
				- i * rect_len(1) /
				frame_points;
			// 1,0;0,1
			rect_frame_V.row(i + 3 * frame_points)
				<< ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
			// 0,0;0,1
		}
		frame_ids = Eigen::VectorXi::LinSpaced(rect_frame_V.rows(), mv_num,
			mv_num + rect_frame_V.rows());
	}

	// Concatenate Vert and Edge
	MatrixXd V;
	MatrixXi E;

	{
		V.resize(V_bnd.rows() + rect_frame_V.rows(), V_bnd.cols());
		V << V_bnd, rect_frame_V;
	}
	E.resize(V.rows(), 2);
	for (int i = 0; i < E.rows(); i++)
		E.row(i) << i, i + 1;
	int acc_bs = 0;
	for (auto bs : bnd_sizes) {
		E(acc_bs + bs - 1, 1) = acc_bs;
		acc_bs += bs;
	}
	E(V.rows() - 1, 1) = acc_bs;
	assert(acc_bs == internal_bnd.size());


	//ycy H 存储的是每个chart的第一个三角形的重心
	MatrixXd H = MatrixXd::Zero(component_sizes.size(), 2);
	{
		int hole_f = 0;
		int hole_i = 0;
		for (auto cs : component_sizes) {
			for (int i = 0; i < 3; i++)
				H.row(hole_i) += m_uv.row(m_T(hole_f, i)); // redoing step 2
			hole_f += cs;
			hole_i++;
		}
	}
	H /= 3.;

	MatrixXd uv2;
	triangulate(V, E, H, uv2, s_T);
	auto bnd_n = internal_bnd.size();

	for (auto i = 0; i < s_T.rows(); i++)
	{
		for (auto j = 0; j < s_T.cols(); j++)
		{
			auto &x = s_T(i, j);
			if (x < bnd_n) x = internal_bnd(x);
			else x += m_uv.rows() - bnd_n;
		}
	}


	{
		w_T.resize(m_T.rows() + s_T.rows(), 3);
		w_T << m_T, s_T;
	}


	w_uv.conservativeResize(m_uv.rows() - bnd_n + uv2.rows(), 2);
	w_uv.bottomRows(uv2.rows() - bnd_n) = uv2.bottomRows(-bnd_n + uv2.rows());
	update_scaffold();
	/*VectorXd uv_max = w_uv.colwise().maxCoeff();
	VectorXd uv_min = w_uv.colwise().minCoeff();*/
	//writeObj(w_uv, surface_F, "G:/a.obj");
}

void ScafData::add_new_patch(const Eigen::MatrixXd &V_in,
                             const Eigen::MatrixXi &F_ref,
                             const Eigen::RowVectorXd &center) {

  const Eigen::MatrixXd& V_ref = V_in;
  Eigen::MatrixXd uv_init;
  Eigen::VectorXi bnd;
  Eigen::MatrixXd bnd_uv;
  std::vector<std::vector<int>> all_bnds;
  boundary_loop(F_ref, all_bnds);
  int num_holes = all_bnds.size() - 1;

  std::sort(all_bnds.begin(), all_bnds.end(), [](auto& a, auto&b){
    return a.size() > b.size();
  });

  bnd =  Map<Eigen::VectorXi>(all_bnds[0].data(),
                              all_bnds[0].size());

  map_vertices_to_circle(V_ref, bnd, bnd_uv);
  bnd_uv *= sqrt(1.0/M_PI);
  bnd_uv.rowwise() += center;
  mesh_measure += 1.0;
  std::cout<<"Mesh Measure "<< mesh_measure <<"; number holes "<< num_holes <<std::endl;

  if(num_holes == 0) 
  {
    if (bnd.rows() == V_ref.rows()) 
	{
      std::cout << "All vert on boundary" << std::endl;
      uv_init.resize(V_ref.rows(), 2);
      for (int i = 0; i < bnd.rows(); i++) {
        uv_init.row(bnd(i)) = bnd_uv.row(i);
      }
    } else 
	{
      Tutte(V_ref.rows(), F_ref, bnd, bnd_uv, uv_init);
    }
  } 
  else 
  {
    auto &primary_bnd = bnd;
    // fill holes
    int n_filled_faces = 0;
    int real_F_num = F_ref.rows();
    for (int i = 0; i < num_holes; i++) 
	{
      n_filled_faces += all_bnds[i + 1].size();
    }
    MatrixXi F_filled(n_filled_faces + real_F_num, 3);
    F_filled.topRows(real_F_num) = F_ref;

    int new_vert_id = V_in.rows();
    int new_face_id = real_F_num;

    for (int i = 0; i < num_holes; i++) 
	{
      int cur_bnd_size = all_bnds[i + 1].size();
      auto it = all_bnds[i + 1].begin();
      auto back = all_bnds[i + 1].end() - 1;
      F_filled.row(new_face_id++) << *it, *back, new_vert_id;
      while (it != back) 
	  {
        F_filled.row(new_face_id++)
            << *(it + 1), *(it), new_vert_id;
        it++;
      }
      new_vert_id++;
    }
    assert(new_face_id == F_filled.rows());
    assert(new_vert_id == V_in.rows() + num_holes);

    Tutte(V_ref.rows()+num_holes,F_filled, primary_bnd, bnd_uv, uv_init);
    uv_init.conservativeResize(V_in.rows(), 2);
  }

  //writeObj(uv_init, F_ref, "tutte.obj");

  component_sizes.push_back(F_ref.rows());

  if (mv_num == 0)
  {
	  w_uv =uv_init;
  }
  else
  {
	  MatrixXd m_uv = w_uv.topRows(mv_num);
	  w_uv.resize(m_uv.rows() + uv_init.rows(), 2);
	  w_uv << m_uv, uv_init;
  }

  m_M.conservativeResize(mf_num + F_ref.rows());
  m_M.bottomRows(F_ref.rows()) = M;

  for(auto cur_bnd : all_bnds) 
  {
    internal_bnd.conservativeResize(internal_bnd.size()+ cur_bnd.size());
    internal_bnd.bottomRows(cur_bnd.size()) = Map<ArrayXi>(cur_bnd.data(),cur_bnd.size()) + mv_num;
    bnd_sizes.push_back(cur_bnd.size());
  }

  m_T.conservativeResize(mf_num + F_ref.rows(), 3);
  m_T.bottomRows(F_ref.rows()) = F_ref.array() + mv_num;
  mf_num += F_ref.rows();

  m_V.conservativeResize(mv_num + V_ref.rows(), 3);
  m_V.bottomRows(V_ref.rows()) = V_ref;
  mv_num += V_ref.rows();

  rect_frame_V.resize(0, 0);

  mesh_improve(true);
}

