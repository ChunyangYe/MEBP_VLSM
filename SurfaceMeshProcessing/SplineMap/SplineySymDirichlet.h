#pragma once

#include "BSplineCurve.h"
#include "BSplineSurface.h"

#include"Parafun.h"
#include"ScafData.h"
#include"PolyMap\SharedData.h"
#include"MeshViewer\MeshDefinition.h"


using namespace Eigen;
using namespace std;

class SplineySymDirichlet
{
public:
	
	SplineySymDirichlet();
	~SplineySymDirichlet();

	void parameterization();
	void load();
	//void write_obj();

	vector<Eigen::Vector3d> samplepoints()
	{
		BSplineSurface	*spline_surface = parafun_solver->spline_surface;
		vector<double> samplepointsx = parafun_solver->samplepointsx;
		vector<double> samplepointsy = parafun_solver->samplepointsy;
		vector<double>	assistx = parafun_solver->assistx;
		vector<double>	assisty = parafun_solver->assisty;
		int Sample_N = samplepointsx.size();
		int Assist_N = assistx.size();
		vector<Eigen::Vector3d> points(Sample_N + Assist_N);
		Eigen::Vector3d p;
		for (int i = 0; i < Sample_N; ++i)
		{
			p = spline_surface->getPoint(samplepointsx[i], samplepointsy[i]);
			p[0] = samplepointsx[i];
			p[1] = samplepointsy[i];
			p[2] = 0;
			points[i] = p;
		}
		for (int i = 0; i < Assist_N; ++i)
		{
			p = spline_surface->getPoint(assistx[i], assisty[i]);
			p[0] = assistx[i];
			p[1] = assisty[i];
			p[2] = 0;
			points[Sample_N + i] = p;
		}
		return points;
	}
	
	vector<vector<Eigen::Vector3d>> & controlpoints()
	{
		BSplineSurface	*spline_surface = parafun_solver->spline_surface;
		return spline_surface->ControlPoints();
	}
	BSplineSurface* splinesurface()
	{
		BSplineSurface	*spline_surface = parafun_solver->spline_surface;
		return spline_surface;
	}

protected:
	std::shared_ptr<Parafun> parafun_solver = nullptr;

	double convgence_con_rate = 1e-5;
	int MAX_ITER_NUM = 500;
	
public:
	//ScafData scaf_data;

};