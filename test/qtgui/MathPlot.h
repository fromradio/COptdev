#ifndef MATH_WIDGET_H__
#define MATH_WIDGET_H__

#include <QStackedWidget>
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_text.h>

class MathPlot : public QStackedWidget
{
	Q_OBJECT

	int 			__plot_size;
	
public:

	MathPlot(const QwtText& text=QwtText("Math Plot"),QWidget *parent = 0);

	QSize minimumSizeHint() const;

	int insertOneQwtPlot();

};

#endif