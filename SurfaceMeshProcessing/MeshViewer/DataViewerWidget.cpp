#include "DataViewerWidget.h"

DataViewerWidget::DataViewerWidget(PolyMap* polymap_, QGLFormat& _fmt, QWidget * _parent)
	:QGLViewerWidget(_fmt, _parent),polymap(polymap_)
{
	srand((unsigned int)time(NULL));
	for (size_t i = 0; i < 100; i++)
	{
		float r_ratio = std::rand() / float(RAND_MAX);
		float g_ratio = std::rand() / float(RAND_MAX);
		float b_ratio = std::rand() / float(RAND_MAX);
		color[i][0] = r_ratio;
		color[i][1] = g_ratio;
		color[i][2] = b_ratio;
	}
}

DataViewerWidget::~DataViewerWidget()
{
}

void DataViewerWidget::draw_scene(int drawmode)
{
	switch (drawmode_y)
	{
	case 0://none

		break;
	case 1://parameterization
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.5f, 2.0f);

		glColor3d(1.0, 1.0, 1.0);
		if (distortion_colors.empty() || distortion_colors.size() != polymap->sD.mesh.n_faces())
		{
			std::cerr << "ERROR: DrawColormap() values error." << std::endl;
			return;
		}
		//glEnable(GL_TEXTURE_2D);
		glBegin(GL_TRIANGLES);
		for (const auto& fh : polymap->sD.mesh.faces())
		{
			glNormal3dv(polymap->sD.mesh.normal(fh).data());
			glColor3d(1.0, 1 - distortion_colors[fh.idx()], 1 - distortion_colors[fh.idx()]);
			for (const auto& fvh : polymap->sD.mesh.fv_range(fh))
			{
				glVertex3dv(polymap->sD.mesh.point(fvh).data());
			}
		}
		glEnd();
		glDisable(GL_POLYGON_OFFSET_FILL);

		draw_src_edges();
		draw_rings(true);
		break;
	case 2://ordermap
		glLineWidth(2);
		glColor3f(0.0, 1.0, 0.0);
		for (size_t i = polymap->hybrid_solver->fix_F_N; i < polymap->hybrid_solver->s_F.rows(); i++)
		{
			glBegin(GL_LINE_LOOP);
			for (size_t j = 0; j < polymap->hybrid_solver->s_F.cols(); j++)
			{
				int vid = polymap->hybrid_solver->s_F(i, j);
				glVertex3d(polymap->hybrid_solver->s_V(vid, 0), polymap->hybrid_solver->s_V(vid, 1), 0.0);
			}
			glEnd();
		}

		draw_src_edges();

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.5f, 2.0f);
		glBegin(GL_TRIANGLES);
		for (size_t i = 0; i < polymap->hybrid_solver->s_F.rows(); i++)
		{
			if (polymap->hybrid_solver->flag_order[i] == 1)
				glColor3d(1.0, 0.7, 0.75);
			else if (polymap->hybrid_solver->flag_order[i] == 2)
				glColor3d(1.0, 0.08, 0.58);
			else
				glColor3d(0.5, 0., 0.5);
			for (size_t j = 0; j < polymap->hybrid_solver->s_F.cols(); j++)
			{
				int vid = polymap->hybrid_solver->s_F(i, j);
				glVertex3d(polymap->hybrid_solver->s_V(vid, 0), polymap->hybrid_solver->s_V(vid, 1), 0.0);
			}
		}
		glEnd();
		glDisable(GL_POLYGON_OFFSET_FILL);

		draw_test();

		//glEnable(GL_POLYGON_OFFSET_FILL);
		//glPolygonOffset(1.5f, 2.0f);
		//glColor3d(1.0, 1.0, 1.0);
		//glBegin(GL_TRIANGLES);
		//for (auto& fh : polymap->sD.mesh.faces())
		//{
		//	glNormal3dv(polymap->sD.mesh.normal(fh).data());
		//	int i = polymap->hybrid_solver->fine2coarse_f[fh.idx()] % 100;
		//	glColor3d(color[i][0], color[i][1], color[i][2]);
		//	for (const auto& fvh : polymap->sD.mesh.fv_range(fh))
		//	{
		//		glVertex3dv(polymap->sD.mesh.point(fvh).data());
		//	}
		//}
		//glEnd();
		//glDisable(GL_POLYGON_OFFSET_FILL);

		break;

	case 3://retriangulation
	{
		//glEnable(GL_POLYGON_OFFSET_FILL);
		//glPolygonOffset(1.5f, 2.0f);
		//glColor3d(1.0, 1.0, 1.0);
		//glBegin(GL_TRIANGLES);
		//for (auto& o1 : polymap->hybrid_solver->order1_fs_b)
		//{
		//	auto fh = polymap->sD.mesh.face_handle(o1);
		//	glNormal3dv(polymap->sD.mesh.normal(fh).data());
		//	glColor3d(1.0, 0.0, 0.0);
		//	for (const auto& fvh : polymap->sD.mesh.fv_range(fh))
		//	{
		//		glVertex3dv(polymap->sD.mesh.point(fvh).data());
		//	}
		//}
		//glEnd();
		//glDisable(GL_POLYGON_OFFSET_FILL);

		//draw_src_edges();

		//for (size_t i = 0; i < polymap->hybrid_solver->boundary_loops_part.size(); i++)
		//{
		//	glLineWidth(5);
		//	glColor3f(color[i][0], color[i][1], color[i][2]);
		//	glBegin(GL_LINE_LOOP);
		//	for (auto&var : polymap->hybrid_solver->boundary_loops_part[i])
		//	{
		//		glVertex3dv(polymap->sD.mesh.point(polymap->sD.mesh.vertex_handle(var)).data());
		//	}
		//	glEnd();
		//}

		//glLineWidth(2);
		//glColor3f(0.0, 1.0, 0.0);
		//for (size_t i = polymap->hybrid_solver->fix_F_N; i < polymap->hybrid_solver->s_F.rows(); i++)
		//{
		//	glBegin(GL_LINE_LOOP);
		//	for (size_t j = 0; j < polymap->hybrid_solver->s_F.cols(); j++)
		//	{
		//		int vid = polymap->hybrid_solver->s_F(i, j);
		//		glVertex3d(polymap->hybrid_solver->s_V(vid, 0), polymap->hybrid_solver->s_V(vid, 1), 0.0);
		//	}
		//	glEnd();
		//}

		draw_test();
		draw_rings();


		//glEnable(GL_POLYGON_OFFSET_FILL);
		//glPolygonOffset(1.5f, 2.0f);
		//glColor3d(1.0, 0.0, 0.0);
		////glEnable(GL_TEXTURE_2D);
		//glBegin(GL_TRIANGLES);
		//auto fh = polymap->sD.mesh.face_handle(4654);
		//glNormal3dv(polymap->sD.mesh.normal(fh).data());
		//for (const auto& fvh : polymap->sD.mesh.fv_range(fh))
		//{
		//	glVertex3dv(polymap->sD.mesh.point(fvh).data());
		//}
		//glEnd();
		//glDisable(GL_POLYGON_OFFSET_FILL);

		break;
	}
	case 4://retriangulation + distortion
		draw_after_iteration();
		break;

	default:
		break;
	}

}

void DataViewerWidget::draw_src_edges()
{
	glLineWidth(0.5);
	glColor3f(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : polymap->sD.mesh.edges())
	{
		auto heh = polymap->sD.mesh.halfedge_handle(eh, 0);
		auto vh0 = polymap->sD.mesh.from_vertex_handle(heh);
		auto vh1 = polymap->sD.mesh.to_vertex_handle(heh);
		glNormal3dv(polymap->sD.mesh.normal(vh0).data());
		glVertex3dv(polymap->sD.mesh.point(vh0).data());
		glNormal3dv(polymap->sD.mesh.normal(vh1).data());
		glVertex3dv(polymap->sD.mesh.point(vh1).data());
	}
	glEnd();
}

void DataViewerWidget::draw_test()
{
	glLineWidth(1);
	glColor3f(0.1, 0.1, 0.1);
	for (size_t i = 0; i < polymap->hybrid_solver->s_F.rows(); i++)
	{
		glBegin(GL_LINE_LOOP);
		for (size_t j = 0; j < polymap->hybrid_solver->s_F.cols(); j++)
		{
			int vid = polymap->hybrid_solver->s_F(i, j);
			glVertex3d(polymap->hybrid_solver->s_V(vid, 0), polymap->hybrid_solver->s_V(vid, 1), 0.0);
		}
		glEnd();
	}

	glPointSize(3);
	glBegin(GL_POINTS);
	glColor3d(0., 0., 0.);
	for (auto&vh : polymap->hybrid_solver->s_mesh.vertices())
	{
		glVertex3dv(polymap->hybrid_solver->s_mesh.point(vh).data());
	}
	glColor3d(1., 1., 0.);
	for (auto&heh : polymap->hybrid_solver->s_mesh.halfedges())
	{
		if (polymap->hybrid_solver->s_mesh.data(heh).var_ids.size() == 0)
			continue;
		OpenMesh::Vec3d v0, v1, v;
		v0 = polymap->hybrid_solver->s_mesh.point(polymap->hybrid_solver->s_mesh.from_vertex_handle(heh));
		v1= polymap->hybrid_solver->s_mesh.point(polymap->hybrid_solver->s_mesh.to_vertex_handle(heh));
		if (polymap->hybrid_solver->s_mesh.data(heh).var_ids.size() == 2)
		{
			v = (2.0*v0 + v1) / 3.0;
			glVertex3dv(v.data());
			v = (v0 + 2.0*v1) / 3.0;
			glVertex3dv(v.data());
		}
		else
		{
			v = (v0 + v1) / 2.0;
			glVertex3dv(v.data());
		}
	}
	glColor3d(1., 0., 0.);
	for (auto&fh : polymap->hybrid_solver->s_mesh.faces())
	{
		if (polymap->hybrid_solver->flag_order[fh.idx()] == 3)
		{
			OpenMesh::Vec3d v(0.,0.,0.);
			for (auto&fvh : polymap->hybrid_solver->s_mesh.fv_range(fh))
			{
				v += polymap->hybrid_solver->s_mesh.point(fvh);
			}
			v /= 3.0;
			glVertex3dv(v.data());
		}
	}
	glEnd();

}

void DataViewerWidget::draw_after_iteration()
{
	//edges
	glLineWidth(1);
	glColor3f(0.0, 1.0, 0.0);
	int N_every = 10;
	double interal_ = 1.0 / N_every;
	for (auto&fh : polymap->hybrid_solver->s_mesh.faces())
	{
		if (fh.idx()>=polymap->hybrid_solver->fix_F_N&&fh.idx()<polymap->hybrid_solver->order1_fix_F_N&&polymap->hybrid_solver->s_mesh.data(fh).var_ids_bases.size()<3)
		{
			glBegin(GL_LINE_LOOP);
			for(const auto&fvh:polymap->hybrid_solver->s_mesh.fv_range(fh))
			{
				glVertex3dv(polymap->hybrid_solver->s_mesh.point(fvh).data());
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_LOOP);
			for (int i = 0; i < N_every; i++)
			{
				double x_ = i * interal_;
				Eigen::Vector3d kesi(1.0 - x_, x_, 0.0);
				auto pos_=polymap->hybrid_solver->calc_img(kesi, polymap->hybrid_solver->s_mesh.data(fh).var_ids_bases);
				glVertex3dv(pos_.data());
			}
			for (int i = 0; i < N_every; i++)
			{
				double x_ = i * interal_;
				Eigen::Vector3d kesi(0.0, 1.0 - x_, x_);
				auto pos_ = polymap->hybrid_solver->calc_img(kesi, polymap->hybrid_solver->s_mesh.data(fh).var_ids_bases);
				glVertex3dv(pos_.data());
			}
			for (int i = 0; i < N_every; i++)
			{
				double x_ = i * interal_;
				Eigen::Vector3d kesi(x_, 0.0, 1.0 - x_);
				auto pos_ = polymap->hybrid_solver->calc_img(kesi, polymap->hybrid_solver->s_mesh.data(fh).var_ids_bases);
				glVertex3dv(pos_.data());
			}
			glEnd();
		}
	}


	////show control point
	//glPointSize(10);
	//glBegin(GL_POINTS);
	//glColor3d(0., 0., 0.);
	////for (int i = 0; i < polymap->hybrid_solver->VAR_N; i++)
	////{
	////	glVertex3d(polymap->hybrid_solver->variables_[i], polymap->hybrid_solver->variables_[i + polymap->hybrid_solver->VAR_N], 0.0);
	////}
	//for (int i = 0; i < polymap->hybrid_solver->s_V.rows()-polymap->sD.frame_ids.size(); i++)
	//{
	//	glVertex3d(polymap->hybrid_solver->variables_[i], polymap->hybrid_solver->variables_[i + polymap->hybrid_solver->VAR_N], 0.0);
	//}
	//glColor3d(1., 1., 0.);
	//for (auto&heh : polymap->hybrid_solver->s_mesh.halfedges())
	//{
	//	if (polymap->hybrid_solver->s_mesh.data(heh).var_ids.size() == 0)
	//		continue;
	//	for (auto&var : polymap->hybrid_solver->s_mesh.data(heh).var_ids)
	//	{
	//		glVertex3d(polymap->hybrid_solver->variables_[var], polymap->hybrid_solver->variables_[var + polymap->hybrid_solver->VAR_N], 0.0);
	//	}
	//}
	//glColor3d(1., 0., 0.);
	//for (auto&fh : polymap->hybrid_solver->s_mesh.faces())
	//{
	//	if (polymap->hybrid_solver->flag_order[fh.idx()] == 3)
	//	{
	//		int var_id = polymap->hybrid_solver->s_mesh.data(fh).var_ids_bases.back().first;
	//		glVertex3d(polymap->hybrid_solver->variables_[var_id], polymap->hybrid_solver->variables_[var_id + polymap->hybrid_solver->VAR_N], 0.0);
	//	}
	//}
	//glEnd();


}

void DataViewerWidget::draw_rings(bool is_afteropt)
{
	if (is_afteropt)
	{
		glLineWidth(3);
		glColor3f(1.0, 1.0, 0.0);
		for (const auto&var : polymap->sD.boundary_loops)
		{
			glBegin(GL_LINE_LOOP);
			for (auto&vid : var)
			{
				glVertex3dv(polymap->sD.mesh.point(polymap->sD.mesh.vertex_handle(vid)).data());
			}
			glEnd();
		}

		int cccc_ = 0;
		for (const auto&var : polymap->hybrid_solver->boundary_loops_part)
		{
			if (std::get<1>(var))
				glColor3f(0.0, 1.0, 0.0);
			else
				glColor3f(0.0, 0.0, 1.0);
				//glColor3f(0.0, (double)cccc_ / polymap->hybrid_solver->boundary_loops_part.size(), 1.0);

			glBegin(GL_LINE_LOOP);
			for (auto&vid : std::get<0>(var))
			{
				glVertex3dv(polymap->sD.mesh.point(polymap->sD.mesh.vertex_handle(vid)).data());
			}
			glEnd();
			cccc_++;
		}
		glLineWidth(5);
		glColor3f(0.5, 1.0, 0.0);
		glBegin(GL_LINE_LOOP);
		for (auto&vid : polymap->sD.boundary_inner_loop)
		{
			glVertex3dv(polymap->sD.mesh.point(polymap->sD.mesh.vertex_handle(vid)).data());
		}
		glEnd();

		glPointSize(20);
		glBegin(GL_POINTS);
		glColor3d(0., 1., 0.);
		for (auto&vhid : polymap->hybrid_solver->engine_seeds_v)
		{
			glVertex3dv(polymap->sD.mesh.point(polymap->sD.mesh.vertex_handle(vhid)).data());
		}
		glEnd();


	}
	else
	{
		glLineWidth(3);
		glColor3f(1.0, 1.0, 0.0);
		for (const auto&var : polymap->sD.boundary_loops)
		{
			glBegin(GL_LINE_LOOP);
			for (auto&vid : var)
			{
				glVertex3d(polymap->hybrid_solver->s_V(polymap->hybrid_solver->fix_v_src2s[vid], 0), polymap->hybrid_solver->s_V(polymap->hybrid_solver->fix_v_src2s[vid], 1), 0.);
			}
			glEnd();
		}

		int cccc_ = 0;
		for (const auto&var : polymap->hybrid_solver->boundary_loops_part)
		{
			if (std::get<1>(var))
				glColor3f(0.0, 1.0, 0.0);
			else
				glColor3f(0.0, (double)cccc_ / polymap->hybrid_solver->boundary_loops_part.size(), 1.0);

			glBegin(GL_LINE_LOOP);
			for (auto&vid : std::get<0>(var))
			{
				glVertex3d(polymap->hybrid_solver->s_V(polymap->hybrid_solver->fix_v_src2s[vid], 0), polymap->hybrid_solver->s_V(polymap->hybrid_solver->fix_v_src2s[vid], 1), 0.);
			}
			glEnd();
			cccc_++;
		}

		glPointSize(10);
		glBegin(GL_POINTS);
		glColor3d(0., 0., 0.);
		for (auto&vhid : polymap->hybrid_solver->engine_seeds_v)
		{
			glVertex3d(polymap->hybrid_solver->s_V(polymap->hybrid_solver->fix_v_src2s[vhid], 0), polymap->hybrid_solver->s_V(polymap->hybrid_solver->fix_v_src2s[vhid], 1), 0.);
		}
		glEnd();

	}
}

void DataViewerWidget::reset_distortion_colormap()
{
	polymap->calc_distortion(SharedData::DIS_WEIGHT::UNIFORM);
	distortion_colors.resize(polymap->sD.distortion.size());
	double d;
	for (size_t i = 0; i < polymap->sD.distortion.size(); i++)
	{
		d = polymap->sD.distortion[i];
		distortion_colors[i] = (d / 4.0 >= 5 ? 0.999999 : (d / 4.0 - 1.0) / 4.0);
	}
}

void DataViewerWidget::set_dataView_drawmode(int dm_)
{
	drawmode_y = dm_;
	switch (dm_)
	{
	case 1://parameterization
		reset_distortion_colormap();
		break;
	case 2://colormap

		break;
	case 3://retriangulation

		break;
	case 4:

		break;
	default:
		break;
	}

}

