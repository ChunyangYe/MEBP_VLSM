#ifndef MESHPROCESSING_MESHPARAMDIALOG_H
#define MESHPROCESSING_MESHPARAMDIALOG_H

#include <QDialog>
#include <QtGui>
#include <QtWidgets>

class MeshParamDialog : public QDialog
{
	Q_OBJECT
public:
	MeshParamDialog(QWidget* parent=0);
	~MeshParamDialog();

	QSize sizeHint()
	{
		QRect rect = QApplication::desktop()->screenGeometry();
		return QSize( int( rect.width()*0.15), rect.height() );
	}
	double get_maxarea()
	{
		return max_area_text->text().toDouble();
	}

private:
	QTabWidget* tabWidget;

signals:
	void init_signal();
	void run_1_iter_signal();
	void autorun_signal();
	void run_cm_signal();
	void smap_signal();
	void hmap_signal();
	void showings_signal(int);
	void solver_meths_signal(int);
	void max_area_signal();
	void load_2dmesh_signal();

private:
	QWidget* Basic_Operation_And_Information;
	QScrollArea *view_BOI;

	QLabel* leftLabel_BOI;

	QGroupBox* g_pipeline;
	QGroupBox* g_one_iteration;
	QGroupBox* g_show;

	QLabel* retri_max_area_label;
	QLineEdit* max_area_text;

	//pipeline
	QPushButton* initialization;
	QPushButton* run_1_iteration;
	QPushButton* load_2d_mesh;
	QPushButton* autorun;
	QPushButton* run_cm;

	//1iteration
	QPushButton* spline_map;
	QPushButton* hybrid_map;

	//showing
	QComboBox* showings_;
	QComboBox* meths_;
	QLabel* show_mode_;

private:
	void create_Basic_Operation_Information_Widget();

private:
	void initDialog();
	void createWidget();
	void createLayout();

};

#endif