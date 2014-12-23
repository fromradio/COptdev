#include <QtGui/QMainWindow>
#include <QtGui/QApplication>

int main( int argc , char *argv[] )
{
	QApplication a(argc,argv);
	QMainWindow  w;
	w.setWindowTitle("Main Window");
	w.show();
	return a.exec();
}