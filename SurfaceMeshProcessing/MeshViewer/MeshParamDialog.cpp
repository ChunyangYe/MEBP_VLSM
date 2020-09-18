#include "MeshParamDialog.h"
#include <QApplication>
#include <QDesktopWidget>

#include "Common\CommonDefinitions.h"

MeshParamDialog::MeshParamDialog(QWidget* parent /* = 0 */)
	:QDialog(parent)
{
	initDialog();
}

MeshParamDialog::~MeshParamDialog()
{
}

void MeshParamDialog::initDialog()
{
	createWidget();
	createLayout();
}

void MeshParamDialog::createWidget()
{
	create_Basic_Operation_Information_Widget();
}

void MeshParamDialog::createLayout()
{
	tabWidget = new QTabWidget();
	tabWidget->addTab(view_BOI, "PolyMap");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(tabWidget, 0, 0, 1, 1);
	setLayout(layout);
}

void MeshParamDialog::create_Basic_Operation_Information_Widget()
{
	g_pipeline = new QGroupBox("Pipeline");
	initialization = new QPushButton("Initialization");
	run_1_iteration = new QPushButton("Run_1_Iteration");
	meths_ = new QComboBox();
	meths_->addItem("SLIM");
	meths_->addItem("CM");
	retri_max_area_label = new QLabel("Retri-MaxArea:");
	max_area_text = new QLineEdit("0.1");
	load_2d_mesh = new QPushButton("Load-2D-Mesh");
	autorun = new QPushButton("AutoRun");
	run_cm = new QPushButton("Run-SRC");
	QGridLayout* layout1 = new QGridLayout();
	int l1_rows = 0;
	layout1->addWidget(load_2d_mesh, l1_rows++, 0, 1, 2);
	layout1->addWidget(initialization, l1_rows++, 0, 1, 2);
	layout1->addWidget(run_1_iteration, l1_rows, 0, 1, 1);
	layout1->addWidget(meths_, l1_rows++, 1, 1, 1);
	layout1->addWidget(retri_max_area_label, l1_rows, 0, 1, 1);
	layout1->addWidget(max_area_text, l1_rows++, 1, 1, 1);
	layout1->addWidget(autorun, l1_rows++, 0, 1, 2);
	layout1->addWidget(run_cm, l1_rows++, 0, 1, 2);
	g_pipeline->setLayout(layout1);

	g_one_iteration = new QGroupBox("Maps");
	spline_map = new QPushButton("SplineMap");
	hybrid_map = new QPushButton("HybridMap");
	QGridLayout* layout2 = new QGridLayout();
	int l2_rows = 0;
	layout2->addWidget(spline_map, l2_rows++, 0, 1, 2);
	layout2->addWidget(hybrid_map, l2_rows++, 0, 1, 2);
	g_one_iteration->setLayout(layout2);

	g_show = new QGroupBox("Show");
	showings_ = new QComboBox();
	showings_->addItem("None");
	showings_->addItem("Parameterization");
	showings_->addItem("OrderMap");
	showings_->addItem("ReTriangulation");
	showings_->addItem("AfterIteration");

	show_mode_ = new QLabel("ShowMode:");

	QGridLayout* layout3 = new QGridLayout();
	int l3_rows = 0;
	layout3->addWidget(show_mode_, l3_rows, 0, 1, 1);
	layout3->addWidget(showings_, l3_rows++, 1, 1, 1);
	g_show->setLayout(layout3);

	leftLabel_BOI = new QLabel("");

	QGridLayout* mainLayout = new QGridLayout();
	int main_index = 0;
	mainLayout->addWidget(g_pipeline, main_index, 0, l1_rows, 1); main_index += l1_rows;
	mainLayout->addWidget(g_one_iteration, main_index, 0, l2_rows, 1); main_index += l2_rows;
	mainLayout->addWidget(g_show, main_index, 0, l3_rows, 1); main_index += l3_rows;

	mainLayout->addWidget(leftLabel_BOI, main_index, 0, 12, 1);


	Basic_Operation_And_Information = new QWidget();
	Basic_Operation_And_Information->setLayout(mainLayout);

	view_BOI = new QScrollArea;
	view_BOI->setFocusPolicy(Qt::NoFocus);
	view_BOI->setFrameStyle(QFrame::NoFrame);
	view_BOI->setWidget(Basic_Operation_And_Information);
	view_BOI->setWidgetResizable(true);

	connect(initialization, SIGNAL(clicked()), SIGNAL(init_signal()));
	connect(run_1_iteration, SIGNAL(clicked()), SIGNAL(run_1_iter_signal()));
	connect(autorun, SIGNAL(clicked()), SIGNAL(autorun_signal()));
	connect(run_cm, SIGNAL(clicked()), SIGNAL(run_cm_signal()));
	connect(spline_map, SIGNAL(clicked()), SIGNAL(smap_signal()));
	connect(hybrid_map, SIGNAL(clicked()), SIGNAL(hmap_signal()));
	connect(showings_, SIGNAL(currentIndexChanged(int)), SIGNAL(showings_signal(int)));
	connect(meths_, SIGNAL(currentIndexChanged(int)), SIGNAL(solver_meths_signal(int)));
	connect(max_area_text, SIGNAL(editingFinished()), SIGNAL(max_area_signal()));
	connect(load_2d_mesh, SIGNAL(clicked()), SIGNAL(load_2dmesh_signal()));
}
