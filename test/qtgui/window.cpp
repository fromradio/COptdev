#include "window.h"



Window::Window(QWidget *parent)
	:
	QMainWindow(parent)
{
	PlotWidget *plotwidget = new PlotWidget(this);
	setCentralWidget(plotwidget);

	connect(plotwidget,SIGNAL(quitSignal()),this,SLOT(close()));
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
