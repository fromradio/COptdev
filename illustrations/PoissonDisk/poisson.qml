import QtQuick 1.0

Rectangle {
	id: page
	width: 500; height: 200
	color: "lightgray"

	Text{
		id:helloText
		text: "Hello World!"
		y: 30
		anchors.horizontalCenter: page.horizontalCenter
		font.pointSize: 24; font.bold: true
	}

	Rectangle{
		width: 5
		height: 5
		x: 20
		y: 20
		color: "green"
		radius: width*0.5
	}

	Grid {
		id: colorPicker
		x: 4; anchors.bottom: page.bottom; anchors.bottomMargin: 4
		rows: 2; columns: 3; spacing: 3

		Cell { cellColor: "red"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "green"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "blue"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "yellow"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "steelblue"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "black"; onClicked: helloText.color = cellColor }
	}
}