

#include <cassert>
#include <cmath>
#include <string>
#include <string>
#include <sstream>

#include <QApplication>
#include <QLabel>
#include <QVBoxLayout>

#include "window.h"


// #include "qwt/qwt_plot.h"
// #include "qwt/qwt_plot_curve.h"
// #include "qwt/qwt_text.h"

#if QWT_VERSION >= 0x060100 || !WIN32
#include "qwt/qwt_point_data.h"
#endif

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);

  Window* window = new Window;
  window->show();

  // QwtPlotCurve * const m_curve = new QwtPlotCurve("Sine");
  // QwtPlot * const m_plot = new QwtPlot(QwtText("CppQwtExample1"));
  // // QwtPlotCurve * const m_

  // m_plot->setStyleSheet("background-color:red;");

  // m_plot->setGeometry(0,0,640,400);
  // m_plot->setAxisScale(QwtPlot::xBottom, 0.0,4.0 * M_PI);
  // m_plot->setAxisScale(QwtPlot::yLeft,-1.5,1.5);
  // std::vector<double> xs;
  // std::vector<double> ys;
  // for (double x = 0; x < 2.0 * M_PI; x+=(M_PI / 50.0))
  // {
  //   xs.push_back(x);
  //   ys.push_back(std::cos(x));
  // }
  // m_curve->setData(new QwtPointArrayData(&xs[0],&ys[0],ys.size()));
  // m_curve->attach(m_plot);
  // m_plot->replot();
  // m_plot->show();
  return a.exec();
}


// #include <QtGui/QMainWindow>
// #include <QtGui/QApplication>
// #include <QtGui/QWidget>
// #include <QtGui/QPushButton>
// #include "window.h"

// int main( int argc , char *argv[] )
// {
// 	QApplication app(argc,argv);
	
	
// 	QWidget widget;
// 	widget.resize(200,120);
// 	QPushButton quit("Quit",&widget);
// 	// quit.resize(75,40);
// 	quit.setFont(QFont("Times",18,QFont::Bold));
// 	QObject::connect(&quit,SIGNAL(clicked()),&app,SLOT(quit()));
// 	quit.setGeometry(10,40,180,40);
// 	QMainWindow  w;
// 	MyPlot	 	p;
// 	w.setWindowTitle("Main Window");
// 	w.show();
// 	p.show();
// 	widget.show();
// 	return app.exec();
// }