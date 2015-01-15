#include "gui.h"
#include "PoissonDisk.h"

PoissonWidget::PoissonWidget(QWidget *parent)
	:
	QWidget(parent)
{
	PoissonDisk pd(2,0.05,10);
	Vector minvec(2);
	minvec[0]=-1;minvec[1]=-1;
	Vector maxvec(2);
	maxvec[0]=1;maxvec[1]=1;
	pd.setRange(minvec,maxvec);
	pd.generate();
	pd.output("test");
}

QSize PoissonWidget::minimumSizeHint() const
{
	return QSize(500,500);
}