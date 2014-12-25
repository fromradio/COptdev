#include "MathPlot.h"

MathPlot::MathPlot(const QwtText& text,QWidget *parent)
	:
	QwtPlot(text,parent)
{
	this->setLineWidth(3);
	this->setFrameStyle(QFrame::Panel|QFrame::Sunken);
}

QSize MathPlot::minimumSizeHint() const
{
	return QSize(500,500);
}