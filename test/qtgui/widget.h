#include <QtGui/QWidget>
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_text.h>


class PlotWidget: public QWidget
{
	Q_OBJECT
public:
	PlotWidget(QWidget *parent = 0 , const QwtText& text = QwtText("COPT display"));
	~PlotWidget(){}
	QSize minimumSizeHint() const;
};
