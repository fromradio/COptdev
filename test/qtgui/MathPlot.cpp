#include "MathPlot.h"

MathPlot::MathPlot(const QwtText& text,QWidget *parent)
	:
	QStackedWidget(parent),
	__text(text),
	__parent(parent),
	__number_of_plots(0)
{
	QwtPlot *plot = new QwtPlot(text,this);
	reinterpret_cast<QFrame*>(plot->canvas())->setFrameStyle(QFrame::NoFrame);
	this->addWidget(plot);
	this->setLineWidth(3);
	this->setFrameStyle(QFrame::Panel|QFrame::Sunken);
}

QSize MathPlot::minimumSizeHint() const
{
	return QSize(500,500);
}

void MathPlot::pushOneQwtPlot()
{
	QwtPlot *plot = new QwtPlot(__text,this);
	reinterpret_cast<QFrame*>(plot->canvas())->setFrameStyle(QFrame::NoFrame);
	this->addWidget(plot);
	__plots.push_back(plot);
	this->setCurrentWidget(plot);
}

void MathPlot::popOneQwtPlot()
{
	QwtPlot *plot = __plots.back();
	this->removeWidget(plot);
	__plots.pop_back();
	delete plot;
	this->setCurrentWidget(__plots.back());
}