#include "SplineySymDirichlet.h"

SplineySymDirichlet::SplineySymDirichlet()
{

}

SplineySymDirichlet::~SplineySymDirichlet()
{
	
}

void SplineySymDirichlet::parameterization()
{
	////load();
	//
	//parafun_solver.reset(new Parafun(sharedd.scaf_data));

	parafun_solver->BPE_spline();

	//std::cout << "use slim" << parafun_solver->use_slim << std::endl;
	//std::cout << "use cm" << parafun_solver->use_cm << std::endl;;
	//std::cout << "use grad" << parafun_solver->use_ngrad << std::endl;
	////write_obj();
}
void SplineySymDirichlet::load()
{
	//scaf_data = ScafData();
	//scaf_data.add_new_patch(V, F, Eigen::RowVector2d(0, 0));
}

//void SplineySymDirichlet::write_obj()
//{
//	Eigen::MatrixXd& V_in = scaf_data.w_uv;
//	int v_num = mesh.n_vertices();
//	for (int i = 0; i < v_num; ++i)
//	{
//		auto vertex = mesh.vertex_handle(i);
//		/*Mesh::Point pos(V_in(i, 0) - 0.5, V_in(i, 1) - 0.5, 0.0);
//		mesh.set_point(vertex, pos/ area_same_factor);*/
//		Mesh::Point pos(V_in(i, 0), V_in(i, 1) , 0.0);
//		mesh.set_point(vertex, pos);
//	}
//
//	Eigen::MatrixXi& F_ref = scaf_data.m_T;
//	string outstr = "./a.obj";
//
//	ofstream of_obj(outstr, ios::trunc);
//
//	if (V_in.cols() == 3)
//	{
//		for (size_t vid = 0; vid < scaf_data.mv_num; vid++)
//		{
//			of_obj << "v " << V_in(vid, 0) << " " << V_in(vid, 1) << " " << V_in(vid, 2) << endl;
//		}
//	}
//	else if (V_in.cols() == 2)
//	{
//		for (size_t vid = 0; vid < scaf_data.mv_num; vid++)
//		{
//			of_obj << "v " << V_in(vid, 0) << " " << V_in(vid, 1) << " " << 0.0 << endl;
//		}
//	}
//
//	for (size_t fi = 0; fi < F_ref.rows(); fi++)
//	{
//		of_obj << "f " << F_ref(fi, 0) + 1 << " " << F_ref(fi, 1) + 1 << " " << F_ref(fi, 2) + 1 << endl;
//	}
//	of_obj.close();
//}
