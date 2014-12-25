# makefile for qttest

all: qtgui
qtgui: qtgui.o window.o window.moc.o widget.moc.o widget.o
	g++ -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt qtgui.o window.o window.moc.o widget.o widget.moc.o -o $@
qtgui.o: qttest.cpp window.h
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt qttest.cpp -o qtgui.o
window.o: window.moc.cpp window.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt window.cpp -o window.o
widget.o: widget.cpp widget.h
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt widget.cpp -o widget.o
window.moc.o: window.moc.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt window.moc.cpp -o $@
widget.moc.o: widget.moc.cpp
	g++ -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui -F/usr/local/lib -framework QtCore -framework QtGui -framework Qwt widget.moc.cpp -o $@
%.moc.cpp: %.h
		moc $< -o $@