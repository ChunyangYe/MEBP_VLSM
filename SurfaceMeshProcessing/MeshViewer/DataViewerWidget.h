#pragma once
#include "QGLViewerWidget.h"
#include <QtWidgets>
#include"PolyMap\PolyMap.h"
class DataViewerWidget :
	public QGLViewerWidget
{
	Q_OBJECT
public:
	DataViewerWidget(PolyMap* polymap_, QGLFormat& _fmt, QWidget* _parent=0);
	~DataViewerWidget();

private:
	void draw_scene(int drawmode);
	void draw_src_edges();
	void draw_test();
	void draw_after_iteration();
	void draw_rings(bool is_afteropt=false);

public:
	void reset_distortion_colormap();
	void set_dataView_drawmode(int dm_);

public:
	PolyMap* polymap;
	int drawmode_y;
	vector<double> distortion_colors;
	float color[100][3];
	
};

