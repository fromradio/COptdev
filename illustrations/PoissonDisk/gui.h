#ifndef POISSON_GUI_H__
#define POISSON_GUI_H__

#include <QWidget>
#include "PoissonDisk.h"
#include <QKeyEvent>

class PoissonWidget
	:
	public QWidget
{
	Q_OBJECT
private:

	PoissonDisk __pd;
public:
	PoissonWidget(double r,int k,QWidget *parent = 0);
	~PoissonWidget(){}

	QSize minimumSizeHint() const;
	QSize sizeHint() const;

	
	void keyPressEvent(QKeyEvent *e);
	void paintEvent(QPaintEvent*);

signals:

public slots:
	void generatePoints();
};

#endif