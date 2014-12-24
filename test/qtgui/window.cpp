#include "window.h"
#include "widget.h"


PlotWidget::PlotWidget(QWidget*parent , const QwtText& text )
	// :
	// QwtPlot(text,parent)
{
	QwtPlot *plotwidget = new QwtPlot(text,this);

	QPushButton *okbutton = new QPushButton("OK",this);
	QGridLayout *gridlayout = new QGridLayout(this);
	
	gridlayout->setMargin(15);
	gridlayout->setSpacing(10);
	gridlayout->addWidget(plotwidget,0,0);
	gridlayout->addWidget(okbutton,1,1);

	plotwidget->replot();
}

QSize PlotWidget::minimumSizeHint() const
{
	return QSize(1000,500);
}

Window::Window(QWidget *parent)
	:
	QMainWindow(parent)
{
	PlotWidget *plotwidget = new PlotWidget(this);
	setCentralWidget(plotwidget);
	// plotwidget->replot();

	// QPushButton *okbutton = new QPushButton("OK",this);
	// QDockWidget *dockwidget = new QDockWidget("test",this);
	// dockwidget->setWidget(okbutton);


	// QPushButton *okbutton = new QPushButton("OK",this);
	// QGridLayout *gridlayout = new QGridLayout(this);
	
	// gridlayout->setMargin(15);
	// gridlayout->setSpacing(10);
	// gridlayout->addWidget(plotwidget,0,0);
	// gridlayout->addWidget(okbutton,1,1);
}
