# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\GeoMiner\mainwindow.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(597, 623)
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        self.centralWidget.setObjectName("centralWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralWidget)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setObjectName("gridLayout")
        self.groupBox = QtWidgets.QGroupBox(self.centralWidget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout_3.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_3.setSpacing(6)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.lineEdit_input = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_input.setMinimumSize(QtCore.QSize(0, 20))
        self.lineEdit_input.setMaximumSize(QtCore.QSize(16777215, 20))
        self.lineEdit_input.setReadOnly(True)
        self.lineEdit_input.setObjectName("lineEdit_input")
        self.gridLayout_3.addWidget(self.lineEdit_input, 0, 0, 1, 1)
        self.pushButton_changeInputDataDir = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_changeInputDataDir.setMinimumSize(QtCore.QSize(80, 21))
        self.pushButton_changeInputDataDir.setMaximumSize(QtCore.QSize(80, 21))
        self.pushButton_changeInputDataDir.setObjectName("pushButton_changeInputDataDir")
        self.gridLayout_3.addWidget(self.pushButton_changeInputDataDir, 0, 1, 1, 1)
        self.lineEdit_output = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_output.setMinimumSize(QtCore.QSize(0, 20))
        self.lineEdit_output.setMaximumSize(QtCore.QSize(16777215, 20))
        self.lineEdit_output.setText("")
        self.lineEdit_output.setReadOnly(True)
        self.lineEdit_output.setObjectName("lineEdit_output")
        self.gridLayout_3.addWidget(self.lineEdit_output, 1, 0, 1, 1)
        self.pushButton_changeOutputDataDir = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_changeOutputDataDir.setMinimumSize(QtCore.QSize(80, 21))
        self.pushButton_changeOutputDataDir.setMaximumSize(QtCore.QSize(80, 21))
        self.pushButton_changeOutputDataDir.setObjectName("pushButton_changeOutputDataDir")
        self.gridLayout_3.addWidget(self.pushButton_changeOutputDataDir, 1, 1, 1, 1)
        self.gridLayout.addWidget(self.groupBox, 0, 0, 1, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralWidget)
        self.groupBox_2.setMaximumSize(QtCore.QSize(16777215, 81))
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.comboBox_GEOid = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_GEOid.setObjectName("comboBox_GEOid")
        self.gridLayout_2.addWidget(self.comboBox_GEOid, 0, 0, 1, 4)
        self.label = QtWidgets.QLabel(self.groupBox_2)
        self.label.setMaximumSize(QtCore.QSize(40, 13))
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 1, 0, 1, 1)
        self.comboBox_pval = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_pval.setMaximumSize(QtCore.QSize(47, 22))
        self.comboBox_pval.setObjectName("comboBox_pval")
        self.gridLayout_2.addWidget(self.comboBox_pval, 1, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox_2)
        self.label_2.setMaximumSize(QtCore.QSize(40, 16777215))
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 1, 2, 1, 1)
        self.comboBox_analysis = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_analysis.setObjectName("comboBox_analysis")
        self.gridLayout_2.addWidget(self.comboBox_analysis, 1, 3, 1, 1)
        self.gridLayout.addWidget(self.groupBox_2, 0, 1, 1, 1)
        self.splitter_2 = QtWidgets.QSplitter(self.centralWidget)
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName("splitter_2")
        self.pushButton_l = QtWidgets.QPushButton(self.splitter_2)
        self.pushButton_l.setObjectName("pushButton_l")
        self.pushButton_clear = QtWidgets.QPushButton(self.splitter_2)
        self.pushButton_clear.setMaximumSize(QtCore.QSize(75, 16777215))
        self.pushButton_clear.setObjectName("pushButton_clear")
        self.pushButton_r = QtWidgets.QPushButton(self.splitter_2)
        self.pushButton_r.setObjectName("pushButton_r")
        self.gridLayout.addWidget(self.splitter_2, 2, 0, 1, 2)
        self.splitter = QtWidgets.QSplitter(self.centralWidget)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.progressBar = QtWidgets.QProgressBar(self.splitter)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.lineEdit_status = QtWidgets.QLineEdit(self.splitter)
        self.lineEdit_status.setMinimumSize(QtCore.QSize(40, 0))
        self.lineEdit_status.setMaximumSize(QtCore.QSize(40, 16777215))
        self.lineEdit_status.setReadOnly(True)
        self.lineEdit_status.setObjectName("lineEdit_status")
        self.pushButton_test_export = QtWidgets.QPushButton(self.splitter)
        self.pushButton_test_export.setMinimumSize(QtCore.QSize(120, 0))
        self.pushButton_test_export.setMaximumSize(QtCore.QSize(120, 16777215))
        self.pushButton_test_export.setObjectName("pushButton_test_export")
        self.gridLayout.addWidget(self.splitter, 3, 0, 1, 2)
        self.tableView = QtWidgets.QTableView(self.centralWidget)
        self.tableView.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.tableView.setSortingEnabled(True)
        self.tableView.setObjectName("tableView")
        self.gridLayout.addWidget(self.tableView, 1, 0, 1, 2)
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 597, 21))
        self.menuBar.setObjectName("menuBar")
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtWidgets.QToolBar(MainWindow)
        self.mainToolBar.setObjectName("mainToolBar")
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtWidgets.QStatusBar(MainWindow)
        self.statusBar.setObjectName("statusBar")
        MainWindow.setStatusBar(self.statusBar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        MainWindow.setTabOrder(self.lineEdit_input, self.pushButton_changeInputDataDir)
        MainWindow.setTabOrder(self.pushButton_changeInputDataDir, self.tableView)
        MainWindow.setTabOrder(self.tableView, self.pushButton_l)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "GEO analyser"))
        self.groupBox.setTitle(_translate("MainWindow", "Data directories"))
        self.pushButton_changeInputDataDir.setText(_translate("MainWindow", "Change input"))
        self.pushButton_changeOutputDataDir.setText(_translate("MainWindow", "Change output"))
        self.groupBox_2.setTitle(_translate("MainWindow", "Experiment"))
        self.label.setText(_translate("MainWindow", "P-value"))
        self.label_2.setText(_translate("MainWindow", "Analysis"))
        self.pushButton_l.setText(_translate("MainWindow", "Case"))
        self.pushButton_clear.setText(_translate("MainWindow", "Clear"))
        self.pushButton_r.setText(_translate("MainWindow", "Control"))
        self.pushButton_test_export.setText(_translate("MainWindow", "Analyse and export"))

