# makefile for qttest

all: qtgui
qtgui: qtgui.o window.o window.moc.o widget.moc.o widget.o mathplot.o mathplot.moc.o
	g++ -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt qtgui.o window.o window.moc.o widget.moc.o widget.o mathplot.o mathplot.moc.o -o $@
qtgui.o: qttest.cpp Window.h
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt qttest.cpp -o qtgui.o
window.o: Window.moc.cpp Window.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt Window.cpp -o window.o
widget.o: Widget.cpp widget.h
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt Widget.cpp -o widget.o
mathplot.o: MathPlot.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt MathPlot.cpp -o $@
window.moc.o: Window.moc.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt Window.moc.cpp -o $@
widget.moc.o: Widget.moc.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt Widget.moc.cpp -o $@
mathplot.moc.o: MathPlot.moc.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt MathPlot.moc.cpp -o $@
%.moc.cpp: %.h
		moc $< -o $@