#ifndef MATH_WIDGET_H__
#define MATH_WIDGET_H__

#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_text.h>

class MathPlot : public QwtPlot 
{
	Q_OBJECT

public:
	MathPlot(const QwtText& text=QwtText("Math Plot"),QWidget *parent = 0);

	QSize minimumSizeHint() const;
	
};

#endif