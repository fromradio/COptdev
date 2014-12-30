#ifndef MATH_WIDGET_H__
#define MATH_WIDGET_H__

#include <QStackedWidget>
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_text.h>
#include <list>

class MathPlot : public QStackedWidget
{
	Q_OBJECT

	QwtText 					__text;
	QWidget* 					__parent;
	std::list<QwtPlot*>			__plots;
	int 						__number_of_plots;
	
public:

	MathPlot(const QwtText& text=QwtText("Math Plot"),QWidget *parent = 0);

	QSize minimumSizeHint() const;

	void pushOneQwtPlot();

	void popOneQwtPlot();

};

#endif