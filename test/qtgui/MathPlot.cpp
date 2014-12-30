#include "MathPlot.h"

MathPlot::MathPlot(const QwtText& text,QWidget *parent)
	:
	QStackedWidget(parent)
{
	this->addWidget(new QwtPlot(text,this));
	this->setLineWidth(3);
	this->setFrameStyle(QFrame::Panel|QFrame::Sunken);
}

QSize MathPlot::minimumSizeHint() const
{
	return QSize(500,500);
}

void QSize 