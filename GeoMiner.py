import pandas as pd
import os
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from PyQt5.QtCore import QAbstractTableModel, Qt, QVariant, QThread, pyqtSignal
from PyQt5.QtGui import QColor
from colour import Color
from geominer_mainwindow import Ui_MainWindow
from glob import glob
import os.path
import GEOparse
import traceback
from scipy.stats import ttest_ind
import json
import re


class PandasModel(QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """

    def __init__(self, df, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = df
        self.colour_these_cols = None
        self.cmap = {}
        self.base_colours = [Color("#53f441"), Color("#cb42f4"), Color("#f45e41"), Color("#418ef4")]
        self.light_colours = [Color("#c5ffbf"), Color("#f6dbff"), Color("#f9c7bd"), Color("#b7d6ff")]
        self.make_colour_palette_lookup_thing()
        self.selected_rows = []
        # self.dataChanged.connect(self.view.refresh)

    def make_colour_palette_lookup_thing(self):
        # First, pick the columns that need pallette lookups in the first place. 
        self.colour_these_cols = self._data[self._data.columns.tolist()[2:-2]]
        self.cmap = {}
        
        i = 0
        for column in self.colour_these_cols.columns.tolist():
            vals = [str(v) for v in self.colour_these_cols[column].unique().tolist()]
            colours = list(self.base_colours[i].range_to(self.light_colours[i], len(vals)))
            self.cmap[column] = {j: k for j, k in zip(vals, colours)}  # I fucking love Python sometimes
            # Take care of looping through colours: 
            # if (for whatever reason) we need more than we have, we'll loop back round.
            i += 1
            if i == len(self.base_colours):
                i = 0
    
    def new_data(self, df):
        # Replace existing data and emit the right signal to make stuff happen
        self._data = df
        self.make_colour_palette_lookup_thing()
        # self.dataChanged.emit()
    
    def rowCount(self, parent=None):
        return len(self._data.values)
    
    def columnCount(self, parent=None):
        return self._data.columns.size
    
    def make_qcolor(self, colour):
        return QColor.fromRgb(*[i*255 for i in colour.rgb])
    
    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            col = index.column()
            row = index.row()
            val = str(self._data.values[row][col])
            if role == Qt.DisplayRole:
                return val
            if role == Qt.BackgroundColorRole:
                # May need to look these up by name?
                colname = self._data.columns[col]
                if colname in self.colour_these_cols.columns.tolist():
                    bgColor = self.make_qcolor(self.cmap[colname][val])
                    return QVariant(QColor(bgColor))
                elif colname == 'sample' and row in self.selected_rows:
                    return QVariant(QColor(self.make_qcolor(Color("green"))))
                else:
                    return QVariant(QColor(self.make_qcolor(Color("white"))))
        return None
    
    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._data.columns[col]
        return None
    
    def set_selected_rows(self, selected_rows):
        self.selected_rows = selected_rows
    
    def set_samples_as_case(self, row_indexes):
        self.layoutAboutToBeChanged.emit()
        for i in row_indexes:
            # Remove mark as control, if any
            self._data.iloc[i, -1] = ""
            # Mark as case
            self._data.iloc[i, 0] = "<"
        self.layoutChanged.emit()
    
    def set_samples_as_control(self, row_indexes):
        self.layoutAboutToBeChanged.emit()
        for i in row_indexes:
            # Remove mark as control, if any
            self._data.iloc[i, 0] = ""
            # Mark as case
            self._data.iloc[i, -1] = ">"
        self.layoutChanged.emit()
    
    def reset_samples_as_case_control(self):
        self.layoutAboutToBeChanged.emit()
        self._data['case'] = ""
        self._data['control'] = ""
        self.layoutChanged.emit()
    
    def get_case_control_samples(self):
        case = self._data[self._data['case'] == "<"].iloc[:, 1].tolist()
        control = self._data[self._data['control'] == ">"].iloc[:, 1].tolist()
        return case, control
    
    def sort(self, Ncol, order):
        self.layoutAboutToBeChanged.emit()
        if order == Qt.DescendingOrder:
            self._data = self._data.sort_values(self._data.columns[Ncol], ascending=False)
        else:
            self._data = self._data.sort_values(self._data.columns[Ncol])
        self.layoutChanged.emit()
    
    def get_rows_with_matching_value(self, index):
        # Ugly as all hell, but it works.
        col = index.column()
        row = index.row()
        val = str(self._data.values[row][col])
        selected_rows = [row]
        
        def get_block_oneway(rowids):
            selected_rows = []
            changed = False
            for i in rowids:
                if self._data.values[i][col] == val and not changed:
                    selected_rows.append(i)
                if self._data.values[i][col] != val:
                    changed = True
            return selected_rows
        
        selected_rows = selected_rows + get_block_oneway(list(reversed(range(0, row))))
        selected_rows = selected_rows + get_block_oneway(list(range(row + 1, len(self._data))))
        return sorted(selected_rows)


class doTTestThread(QThread):
    progbar = pyqtSignal(int)
    
    def __init__(self, sample_name, case, control, threshold_pval, analysis, data, settings):
        QThread.__init__(self)
        self.sample_name = sample_name
        self.case = case
        self.control = control
        self.threshold_pval = threshold_pval
        self.analysis = analysis
        self.data = data
        self.settings = settings

    def run(self):
        try:
            genes_up = []
            genes_down = []
            if self.analysis in [1, 2]:
                self.threshold_pval = self.threshold_pval / 2
            c = 0
            for i, row in self.data.iterrows():
                c += 1
                if c % 100 == 0:
                    self.progbar.emit(int((c / len(self.data)) * 100))
                gene = str(row['IDENTIFIER'])
                casevals = [float(row[d]) for d in self.case if row[d] != 'null']
                controlvals = [float(row[d]) for d in self.control if row[d] != 'null']
                t, prob = ttest_ind(casevals, controlvals)
                if prob <= self.threshold_pval:
                    # print(gene + "\t" + str(t) + "\t" + str(prob))
                    if t > 0 and self.analysis in [0, 1]:
                        genes_up.append(gene)
                    elif t < 0 and self.analysis in [0, 2]:
                        genes_down.append(gene)
            self.progbar.emit(100)
            # Now let's write our output files.
            # outdir = self.settings['outdir'] + "/" + re.sub(".soft", "", self.sample_name)
            # if not os.path.exists(outdir):
            #     os.mkdir(outdir)
            # with open(outdir + "/up.txt", 'w+') as file_handler:
            #     for item in genes_up:
            #         file_handler.write("{}\n".format(item))
            # with open(outdir + "/down.txt", 'w+') as file_handler:
            #     for item in genes_down:
            #         file_handler.write("{}\n".format(item))

            # Actually, let's just have the one output file. 
            frame_up = pd.DataFrame([[i, 1.0] for i in genes_up], columns=['gene', 'weight'])
            frame_down = pd.DataFrame([[i, -1.0] for i in genes_down], columns=['gene', 'weight'])
            frame_out = pd.concat([frame_up, frame_down])
            frame_out.to_csv(self.settings['outdir'] + "/" + re.sub(".soft", ".csv", self.sample_name), index=False)
        except Exception as e:
            print("ERROR IN TEST THREAD")
            traceback.print_exc()
        

class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)
        self.make_ready()
        self.softfiles = []
        self.test_thread = None
        self.fac = None
        self.data = None
        self.selected_experimental_level = []
        
        self.pushButton_changeInputDataDir.clicked.connect(self.open_directory_of_geo_files)
        self.pushButton_changeOutputDataDir.clicked.connect(self.open_output_directory)
        self.comboBox_GEOid.currentIndexChanged.connect(self.get_experimental_design)
        
        self.pushButton_l.clicked.connect(self.set_case)
        self.pushButton_r.clicked.connect(self.set_control)
        self.pushButton_test_export.clicked.connect(self.do_test)
        
        self.comboBox_pval.addItems(['0.10', '0.05', '0.01'])
        self.comboBox_pval.setCurrentIndex(2)
        
        self.comboBox_analysis.addItems(['A <> B', 'A > B', 'A < B'])
        
        self.show()
        self.settings = {}
        self.boot_up()
    
    def boot_up(self):
        with open('GeoMiner_settings.json', 'r') as file:
            self.settings = json.load(file)
        self.lineEdit_input.setText(self.settings['datadir'])
        self.lineEdit_output.setText(self.settings['outdir'])
        self.read_directory_of_geo_files()
    
    def make_ready(self):
        self.progressBar.setValue(0)
        self.lineEdit_status.setText("Ready")
    
    def open_output_directory(self):
        self.settings['outdir'] = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.lineEdit_output.setText(self.settings['outdir'])
        with open('GeoMiner_settings.json', 'w') as outfile:
            json.dump(self.settings, outfile)
    
    def open_directory_of_geo_files(self):
        self.settings['datadir'] = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.lineEdit_input.setText(self.settings['datadir'])
        with open('GeoMiner_settings.json', 'w') as outfile:
            json.dump(self.settings, outfile)
        self.read_directory_of_geo_files()
    
    def read_directory_of_geo_files(self):
        self.softfiles = [os.path.basename(i) for i in glob(self.settings['datadir'] + "/*.soft")]
        self.comboBox_GEOid.addItems(self.softfiles)
    
    def get_experimental_design(self):
        self.make_ready()
        experiment_id = self.comboBox_GEOid.currentText()
        print(experiment_id)
        # OK, this is where things start getting complicated. First, we've got to read the .soft file, and then
        # we've got to fit it into a table model.
        gse = GEOparse.get_GEO(filepath=self.settings['datadir'] + "/" + experiment_id)
        print(gse)
        # Looks like we want gse.columns. How can I make that a dataframe?
        # Oh - it already is one. Boss. 
        self.fac = gse.columns
        # gsm.table is the actual values, conveniently already as a pandas dataframe.  
        self.data = gse.table
        # Boss. Now we can set up a table view and all that good shit.
        self.fac['sample'] = self.fac.index
        self.fac['case'] = ''
        self.fac['control'] = ''
        mdl = PandasModel(self.fac[['case', 'sample'] +
                                   [i for i in self.fac.columns if
                                    i not in ['sample', 'case', 'control', 'description']] +
                                   ['sample', 'control']])
        self.tableView.setModel(mdl)
        self.tableView.resizeColumnsToContents()
        # This has got to be done here, it seems, because it won't work until after a model is set. 
        self.tableView.selectionModel().selectionChanged.connect(self.deal_with_selecting_groups)
        self.pushButton_clear.clicked.connect(self.tableView.model().reset_samples_as_case_control)
    
    def deal_with_selecting_groups(self):
        self.make_ready()
        # When you click on one thing, we should go and get indices up and down the column until we run into something 
        # different. Then we should auto-select them all. 
        inds = self.tableView.selectionModel().selectedIndexes()
        if inds:
            blox = []
            for i in inds:
                blox.extend(self.tableView.model().get_rows_with_matching_value(i))
            self.selected_experimental_level = sorted(list(set(blox)))
            # Fuckin A. We'll leave it at that - otherwise this quickly gets very complicated.
    
    def set_case(self):
        self.make_ready()
        if self.selected_experimental_level:
            self.tableView.model().set_samples_as_case(self.selected_experimental_level)
            self.selected_experimental_level = []
            self.tableView.clearSelection()
    
    def set_control(self):
        self.make_ready()
        if self.selected_experimental_level:
            self.tableView.model().set_samples_as_control(self.selected_experimental_level)
            self.selected_experimental_level = []
            self.tableView.clearSelection()
    
    def update_progbar(self, i):
        self.progressBar.setValue(i)
    
    def test_done(self):
        self.lineEdit_status.setText("Done!")
    
    def do_test(self):
        # Alright, now let's get down to business. 
        # First: get sample names marked as case and control.
        case, control = self.tableView.model().get_case_control_samples()
        self.test_thread = doTTestThread(sample_name=self.comboBox_GEOid.currentText(),
                                         case=case, 
                                         control=control,
                                         threshold_pval=float(self.comboBox_pval.currentText()),
                                         analysis=self.comboBox_analysis.currentIndex(),
                                         data=self.data,
                                         settings=self.settings)
        self.test_thread.progbar.connect(self.update_progbar)
        self.test_thread.finished.connect(self.test_done)
        self.lineEdit_status.setText("Testing")
        self.test_thread.start()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWin = MainWindow()
    ret = app.exec_()
    sys.exit(ret)
