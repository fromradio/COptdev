#include "widget.h"
#include <QVBoxLayout>
#include "qwt/qwt_point_data.h"
#include <iostream>
#include <algorithm>
#include <QFileDialog>
#include <fstream>
#include <QPalette>
#include <QColorDialog>


PlotWidget::PlotWidget(QWidget*parent , const QwtText& text )
	:
	__plot_widget(new MathPlot(text,this)),
	__color_dialog(new QColorDialog(this))
{
	// __plot_widget->canvas()->setStyleSheet(QString::fromUtf8("border:1px solid red"));
	__color_dialog->hide();
	__color_dialog->setCurrentColor(__plot_widget->currentWidget()->palette().color(QPalette::Background));
	// __plot_widget->setMargin(10);
	// reinterpret_cast<QFrame*>(__plot_widget->currentWidget()->canvas())->setFrameStyle(QFrame::NoFrame);
	QGridLayout *gridlayout = new QGridLayout(this);
	
	gridlayout->setMargin(30);
	gridlayout->setSpacing(10);
	gridlayout->addWidget(__plot_widget,0,0);

	QVBoxLayout *vlayout 		= new QVBoxLayout;
	QPushButton *testbutton 	= new QPushButton("Test");
	QPushButton *readbutton 	= new QPushButton("Read Data");
	QPushButton *cancelbutton 	= new QPushButton("Cancel");
	QPushButton *changebutton 	= new QPushButton("Background");

	/** set the vlayout */
	vlayout->addStretch();
	vlayout->addWidget(testbutton);
	vlayout->addWidget(readbutton);
	vlayout->addWidget(changebutton);
	vlayout->addWidget(cancelbutton);

	gridlayout->addLayout(vlayout,0,1);

	connect(cancelbutton,SIGNAL(clicked()),this,SLOT(quitSlot()));
	connect(testbutton,SIGNAL(clicked()),this,SLOT(testPlot()));
	connect(readbutton,SIGNAL(clicked()),this,SLOT(openFile()));
	connect(changebutton,SIGNAL(clicked()),this,SLOT(changeBackground()));

	connect(__color_dialog,SIGNAL(currentColorChanged(const QColor&)),this,SLOT(setPlotWidgetBackground(const QColor&)));

	reinterpret_cast<QwtPlot*>(__plot_widget->currentWidget())->replot();
}

QSize PlotWidget::minimumSizeHint() const
{
	return QSize(500,500);
}

void PlotWidget::quitSlot()
{
	emit(quitSignal());
}

void PlotWidget::changeBackground()
{
	if (!__color_dialog->isVisible())
		__color_dialog->show();
}

void PlotWidget::setPlotWidgetBackground(const QColor& col)
{
	QPalette pal(__plot_widget->palette());
	pal.setColor(QPalette::Background,col);
	__plot_widget->currentWidget()->setAutoFillBackground(true);
	__plot_widget->currentWidget()->setPalette(pal);
	reinterpret_cast<QwtPlot*>(__plot_widget->currentWidget())->replot();
}

void PlotWidget::openFile()
{
	QString filename = QFileDialog::getOpenFileName(this,tr("Open files"),"./",tr("COPT Files(*.copt)"));
	if (filename.isEmpty())
		return;
	std::ifstream fin(filename.toStdString().c_str());
	std::vector<double> vals;
	int size;
	fin>>size;
	vals.reserve(size);
	do{
		double temp;
		if(fin>>temp)
			vals.push_back(temp);
	}while(fin);
	std::vector<double> indices(vals.size());
	std::cout<<vals.size()<<std::endl;
	for ( int i = 0 ; i < vals.size() ; ++ i )
		indices[i] = i;
	readData(&indices[0],&vals[0],vals.size());
}

void PlotWidget::readData( const double *xdata, const double *ydata, int size )
{
	QwtPlotCurve *const curve = new QwtPlotCurve;
	curve->setSamples(xdata,ydata,size);
	curve->attach(reinterpret_cast<QwtPlot*>(__plot_widget->currentWidget()));
	curve->setPen(QColor(255,0,0),3.0);
	double minx = *std::min_element(xdata,xdata+size);
	double maxx = *std::max_element(xdata,xdata+size);
	double miny = *std::min_element(ydata,ydata+size);
	double maxy = *std::max_element(ydata,ydata+size);
	reinterpret_cast<QwtPlot*>(__plot_widget->currentWidget())->setAxisScale(QwtPlot::xBottom,minx,maxx);
	reinterpret_cast<QwtPlot*>(__plot_widget->currentWidget())->setAxisScale(QwtPlot::yLeft,miny,maxy);
	reinterpret_cast<QwtPlot*>(__plot_widget->currentWidget())->replot();

}

void PlotWidget::testPlot()
{
	// std::vector<double> xs;
	// std::vector<double> ys;
	// QwtPlotCurve *const curve = new QwtPlotCurve("Sine");
	// int n = 100;
	// xs.resize(n);
	// ys.resize(n);
	// for ( int i = 0 ; i < n ; ++ i )
	// {
	// 	xs[i] = static_cast<double>(1.0)/(n-1)*i;
	// 	ys[i] = xs[i]*xs[i];
	// }
	// this->readData(&xs[0],&ys[0],n);
	setPlotWidgetBackground(QColor(0,255,0));
	// curve->setPen(QColor(255,0,0),3.0);
	// curve->setSamples(&xs[0],&ys[0],n);
	// curve->attach(__plot_widget);
	// __plot_widget->setAxisScale(QwtPlot::xTop,0.0,1.0);
	// __plot_widget->setAxisScale(QwtPlot::yLeft,0.0,1.0);
	// __plot_widget->replot();
}