#ifndef WINDOW_H__
#define WINDOW_H__

#include <QtGui/QWidget>
#include <QtGui/QMainWindow>
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_text.h>
#include <QtGui/QPushButton>
#include <QtGui/QGridLayout>
#include <QDockWidget>

class DrawWidget : public QWidget
{
	Q_OBJECT
public:
	DrawWidget(QWidget *parent = 0);
};

class Window
	:
	public QMainWindow
{
	Q_OBJECT
public:
	Window( QWidget *parent = 0 );
	~Window(){}
};




// class MyPlot: public QwtPlot
// {
// public:
// 	MyPlot(char *name = 0 ,QWidget *parent = 0):QwtPlot(name,parent)
// 	{
// 		// Show a title
// 		setTitle("This is an Example ");

// 		// Show a legend at the bottom
// 		// setAutoLegend(true);
// 		setLegendPos(Qwt::Bottom);

// 		// Show the axes
// 		setAxisTitle(xBottom,"x");
// 		setAxisTitle(yLeft, "y");

// 		// Insert two curves and get IDs for them
// 		long cSin = insertCurve("y=sin(x)");
// 		long cSign = insertCurve("y=sign(sin(x))");

// 		// Calculate the data, 500 points each
// 		const int points = 500;
// 		double x[points];
// 		double sn[points];
// 		double sg[points];

// 		for ( int i = 0 ; i < points ; ++ i )
// 		{
// 			x[i] = (3.0*3.14/double(points))*double(i);

// 			sn[i] = 2.0*sin(x[i]);
// 			sg[i] = (sn[i]>0)?1:((sn[i]<0)?-1:0);
// 		}

// 		// Copy the data to the plot
// 		setCurveData( cSin, x, sn, points );
// 		setCurveData( cSign, x, sg, points );

// 		// set the style of the curves
// 		setCurvePen( cSin, QPen(Qt::blue));
// 		setCurvePen( cSign, QPen(Qt::green,3));

// 		// show the plots
// 		replot();
// 	}
// };


#endif
