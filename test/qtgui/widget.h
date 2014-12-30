#ifndef WIDGET_H__
#define WIDGET_H__

#include <QtGui/QWidget>
// #include <qwt/qwt_plot.h>
// #include <qwt/qwt_plot_curve.h>
// #include <qwt/qwt_text.h>
#include <QtGui/QPushButton>
#include <QtGui/QGridLayout>
#include <QColorDialog>
#include "MathPlot.h"


class PlotWidget: public QWidget
{
	Q_OBJECT

	/** the plot widget */
	MathPlot			*const __plot_widget;

	QColorDialog 		* __color_dialog;
public:
	PlotWidget(QWidget *parent = 0 , const QwtText& text = QwtText("COPT display"));
	~PlotWidget(){}
	QSize minimumSizeHint() const;

signals:
	void quitSignal();
public slots:
	void quitSlot();

	/** open a COPT data file */
	void openFile();

	void readData( const double *xdata, const double *ydata, int size );

	/** for debug */
	void testPlot();

	/** change background color */
	void changeBackground();

	/** set the background color of plot widget */
	void setPlotWidgetBackground(const QColor& col);
};

#endif