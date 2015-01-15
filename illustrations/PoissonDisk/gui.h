#ifndef POISSON_GUI_H__
#define POISSON_GUI_H__

#include <QWidget>

class PoissonWidget
	:
	public QWidget
{
	Q_OBJECT
private:

public:
	PoissonWidget(QWidget *parent = 0);
	~PoissonWidget(){}

	QSize minimumSizeHint() const;

signals:
public slots:
};

#endif