#Import libraries
import pandas as pd
import os
import tkinter
import time
import functools
from tkinter import filedialog
import sys
from PyQt5.QtWidgets import (QFileDialog, QAbstractItemView, QListView,
                            QTreeView, QDialog, QApplication)
from PyQt5.QtCore import Qt as QtCore_Qt
from PyQt5 import Qt as PyQt5_Qt
import datetime
import pickle
import itertools

"""
FUNCTIONS
"""

#It's important that all the pickles used are saved using the pickle_save function, to allow the reading with the
#pickle_load function. Also it detects incompatibilities with the versions

def printnumberedlist(list):
    """
    Takes all the elements in a list and puts them in a single
    string, per line, and appending a number before each one.

    The  string is meant to be printed, as a way to show the elements
    of the list.

    Parameters
    ----------
    list : list
        A non nested list.

    Returns
    -------
    str
        A string containing all the list elements, prepared to be printed.

    """
    counter = 0
    output_str = ''
    for e in list:
        output_str = ''.join((output_str, '(', str(counter), ') ',
                              str(e), os.linesep))
        counter += 1
    return output_str

def listupdater(new, old):
    """
    Updates a list (old) with the elements of another one (new).
    The resulting one won't have repeated elements.

    Parameters
    ----------
    new : list
        A list with the elements to be added.
    old : list
        The list that will be updated with the elements of the other list.

    Returns
    -------
    list
        A sorted list combining both lists, without repeated elements.

    """
    if new != None:
        list_set = set(old)
        list_set.update(new)
        old= sorted(list(list_set))
    return old


def fileSelector(title='Select files', firstdir='', filetypes=''):
    """
    Select one or more files.

    Parameters
    ----------
    title : str, optional
        A string that will appear as the title of the window.
        Default is: 'Select files'.
    firstdir : str, optional
        The first path that will be explored. Default is: '', the
        directory where the script is located or the one where the last file
        was selected with this function.
    filetypes : array_like, optional
        A vector containing description-extension pairs. Leaving this parameter
        empty "()" or with a "''" disables this feature. The script expects at
        least two file types. The extension filtering is disabled by default.

        To add a two-types filter, an argument like the following must be
        entered:
            filetypes=(('Jpeg files', '*.jpg'), ('All files', '*.*'))

        To add a file type with two or more possible extensions:
            filetypes=(('Jpeg files', ('*.jpg', '*.jpeg')),
            ('All files', '*.*'))
        Or:
            filetypes=(('Jpeg files', '*.jpg'), ('Jpeg files', '*.jpeg'))

        The function can be tricked to show only one required file type and
        extension:
            filetypes=(('Jpeg files', '*.jpg'), ('Jpeg files', ''))

    Returns
    -------
    file_list : list or None
        A list containing the paths to the chosen files is returned if at least
        one file has been selected.
        None is returned if any file has been selected.

    """
    path = tkinter.Tk()
    path.withdraw()
    path.wm_attributes('-topmost', 1)
    file_list = filedialog.askopenfilenames(parent=path, title=title,
                                            initialdir=firstdir, filetypes=filetypes)
    path.destroy()
    if file_list == '':
        return None
    else:
        return list(file_list)


def folderSelector():
    """
    Select one or more folders.

    Parameters
    ----------
    None

    Returns
    -------
    list or None
        A list containing the paths to the chosen folders is returned if at
        least one folder has been selected.
        None is returned if any folder has been selected.

    """
    class getExistingDirectories(QFileDialog):
        def __init__(self, *args):
            super(getExistingDirectories, self).__init__(*args)
            self.setOption(self.DontUseNativeDialog, True)
            self.setFileMode(self.Directory)
            self.setOption(self.ShowDirsOnly, True)
            self.findChildren(QListView)[0].setSelectionMode(
                QAbstractItemView.ExtendedSelection)
            self.findChildren(QTreeView)[0].setSelectionMode(
                QAbstractItemView.ExtendedSelection)
            self.show()
            self.setWindowFlags(QtCore_Qt.WindowStaysOnTopHint)
    qapp = QApplication(sys.argv)
    dlg = getExistingDirectories()
    if dlg.exec_() == QDialog.Accepted:
        qapp.exit()
        return dlg.selectedFiles()


def filefolderManager(mode, fileParams={}, msg=None):
    """
    Handles file and folder selection to retrieve a list of elements.
    Note that it returns the absolute path of each selected item.

    Use msg to add messages and instructions to
    the main menu and to the selection dialog respectively.

    Parameters
    ----------
    mode : str
        The desired target of the function. There's two accepted modes: 'file'
        and 'folder'.
    fileParams : dict, optional
        A dictionary containing the keyword arguments that will be unpacked to
        fileSelector. The dictionary is empty by default (any argument
        will be passed to the function).
    msg : str or None, optional
        An string to be printed before the description of the main menu options
        (and after listing the already chosen items after at least one has been
        selected). A None value disables this feature. The default value is
        None.

    Returns
    -------
    data_dict: dict or None
        A dictionary of the data selected in which each key is a group
        of your experiment, it has to be located in the first element of
        the name of your file,
        sorted list of the chosen elements (files or folders). If any
        file or folder is selected or there's a problem with the mode, it
        returns None.

    """
    elemList = []
    print_dict = {'file': {'menu1': ''.join(('Select an option:',
                                             os.linesep, '(s) Select files.', os.linesep, '(e) Exit.')),
                           'menu2': ''.join(('Select multiple files or click "Cancel" to stop.')),
                           'opt1': ''.join(('You can remove specific files by entering ',
                                            'the number next to its name.'))},
                  'folder': {'menu1': ''.join(('Please select an option:',
                                               os.linesep, '(s) Select folder', os.linesep,
                                               '(e) Exit.')),
                             'menu2': ''.join(('Select folder or click "Cancel" to stop.')),
                             'opt1': ''.join(('You can also remove folder by ',
                                              'entering the number next to its name.'))}}
    def_dict = {'file': functools.partial(fileSelector,
                                          **fileParams), 'folder': folderSelector}
    menuinput = None
    if mode not in def_dict.keys():
        menuinput = 'e'
        print(''.join(('[Error]: Wrong mode, only file or folders mode are accepted')))
    while menuinput != 'e':
        if len(elemList) == 0 and menuinput != 's':
            if not msg == None: print(msg, end='')
            menuinput = input(os.linesep.join((print_dict[mode]['menu1'],
                                               'Option: '))).lower()
        elif menuinput == 's':
            print(print_dict[mode]['menu2'])
            elemList = listupdater(def_dict[mode](), elemList)
            menuinput = None
            if len(elemList) == 0:
                continue
        elif menuinput in list(map(str, range(len(elemList)))):
            menuinput = int(menuinput)
            del elemList[menuinput]
            if len(elemList) == 0:
                menuinput = None
                continue
            print(printnumberedlist(elemList), end='')
            if not msg == None: print(os.linesep + msg, end='')
            menuinput = input(os.linesep.join((print_dict[mode]['menu1'],
                                               print_dict[mode]['opt1'], 'Option: '))).lower()
        elif menuinput not in list(map(str, range(len(elemList)))) and \
                menuinput != 'a':
            print(printnumberedlist(elemList), end='')
            if not msg == None: print(os.linesep + msg, end='')
            menuinput = input(os.linesep.join((print_dict[mode]['menu1'],
                                               print_dict[mode]['opt1'], 'Option: '))).lower()
    if len(elemList) == 0:
        elemList = None
    else:
        elemList = [e.replace('/', os.path.sep) for e in elemList]
    return elemList


def pickle_save(variable, filepath):
    """
    Serializes a variable using pickle.

    Note that the serialization format is guaranteed to be backwards compatible
    across Python releases. Anyways, pickle and Python version number is packed
    together with the chosen variable, to be checked when the file is loaded.

    Parameters
    ----------
    variable : variable
        Initialized variable.
    filepath : str
        Path to the file that will be created. As it uses the pickle
        serialization format it must use the '.pickle' extension.
        E.g. 'my_variable.pickle'.

    Returns
    -------
    None

    """
    python_version = '.'.join([str(sys.version_info[i]) for i in range(3)])
    pickle_version = pickle.format_version
    save_dict = {'version_dict': {'Python': python_version,
                                  'pickle': pickle_version}, 'variable': variable}
    with open(filepath, 'wb') as fileopen:
        pickle.dump(save_dict, fileopen)
    return None


def pickle_load(filepath, version_warning=True):
    """
    Loads a variable using pickle.

    It also optionally outputs a warning if there's a version mismatch between
    packages used, comparing the saving and loading processes.

    Parameters
    ----------
    filepath : str
        Path to the file that will be loaded.
    version_warning : bool, optional
        Boolean to control whether a warning is displayed if there's a version
        mismatch between the packages used when the pickle file was created
        against the ones used when loading the data. Default value is True
        (warning enabled).

    Returns
    -------
    variable : variable
        The element stored in the pickle file.

    """
    with open(filepath, 'rb') as fileopen:
        load_dict = pickle.load(fileopen)
    variable = load_dict['variable']
    version_dict = load_dict['version_dict']
    python_version = '.'.join([str(sys.version_info[i]) for i in range(3)])
    pickle_version = pickle.format_version
    equal_version = version_dict['Python'] == python_version and version_dict[
        'pickle'] == pickle_version
    if version_warning and not equal_version:
        warnings.warn(''.join(('Dependency version mismatch found. Current ',
                               'versions: Python==', python_version, ', pickle==',
                               pickle_version, '. Chosen file versions: Python==',
                               version_dict['Python'], ', pickle==', version_dict['pickle'],
                               '.')), category=Warning)
    return variable

def data_wrangling(filepath, bilateral = False, ROI = False):
    """
    Selection, sorting and data preparation of the csv from lifecanvas
    to create another dataframe or csv file for the network analyses

    Parameters
    ----------
    filepath: string containing the path were the csv files and the ROI file are located
    Bilateral: Boolean that defines if you can to do bilateral analyses or all combine.
                True: separate left and right
                False: mean of left and right
    ROI: Boolean if you want to organize your data anatomically
    Returns
    -------
    dict or None
        A dictionary containing the data_wrangled dataframes per condition
        in your experiment (keys). All the data is saved in yoyr folder as a csv
        None is returned if any file has been selected.

    """
    files = [filepath + '/' + f for f in os.listdir(filepath) if f.endswith('.csv') or f.endswith('.pickle')]
    data_dict = {}
    discard = False
    for i in files:
        if os.path.basename(i).split('.')[0] == 'ROIs':
            ROIs = pd.read_csv(i)
            ROIs = ROIs.loc[:, ["Abbreviation", "Allen Group Name"]].sort_values("Allen Group Name").set_index(
                "Abbreviation").T
            cols = ROIs.columns.tolist()
            if bilateral:
                cols = [['-'.join([c,'L']), '-'.join([c,'R'])] for c in cols]
                cols = list(itertools.chain.from_iterable(cols))
        elif "Discard" in os.path.basename(i).split('.')[0]:
            discard = pickle_load(i)
            print("Selected areas will be discarded")
        elif i.endswith('.pickle'):
            print('The selected folder already contains a data wrangled csv file')
        elif not 'parent_structure_id' in list(pd.read_csv(i).columns):
            print('The selected folder already contains a data wrangled csv file')
        else:
            if os.path.basename(i).split('.')[0].split('_')[0] not in data_dict.keys():
                data_dict[os.path.basename(i).split('.')[0].split('_')[0]] = [pd.read_csv(i).assign(filename=lambda \
                                x: os.path.basename(i).split('.')[0], condition=lambda \
                                x: os.path.basename(i).split('.')[0].split('_')[0],animal = lambda x: \
                                os.path.basename(i).split('.')[0].split('_')[1])]
            else:
                data_dict[os.path.basename(i).split('.')[0].split('_')[0]].append(pd.read_csv(i).assign(filename=lambda \
                                x: os.path.basename(i).split('.')[0], condition=lambda \
                                x: os.path.basename(i).split('.')[0].split('_')[0],animal = lambda x: \
                                os.path.basename(i).split('.')[0].split('_')[1]))
    data = {exp: pd.concat(df, ignore_index = True).loc[:,["name","acronym","density (cells/mm^3)","filename", "condition", "animal"]]\
            for exp,df in data_dict.items()}
    if discard:
        data_discard = {exp:df[~df.name.str.contains('|'.join(discard))] for exp,df in data.items()}
        data_deleted = {exp:df[~df.name.str.islower()] for exp,df in data_discard.items()}
    else:
        data_deleted = {exp:df[~df.name.str.islower()] for exp,df in data.items()}
    data_prewrangled = {}
    for exp,df in data_deleted.items():
        df = df.copy()
        if not bilateral:
            split_bool = df.name.str.contains("left")
            left_df = df[split_bool]
            left_df["name"] = left_df["name"].map(lambda x: x.lstrip('left'))
            right_df = df[~split_bool]
            right_df["name"] = right_df["name"].map(lambda x: x.lstrip('right'))
            merged_df = pd.merge(left_df, right_df, how="inner", on=["name", "filename", "animal", "condition"], suffixes=("_Left", "_Right"))
            merged_df["Bilateral Density (cells/mm^3)"] = merged_df[["density (cells/mm^3)_Left", "density (cells/mm^3)_Right"]].mean(
                            axis=1)
            merged_df = merged_df[["filename","condition", "animal", "name", "acronym_Left", "Bilateral Density (cells/mm^3)"]]
            merged_df["acronym_Left"] = merged_df["acronym_Left"].map(lambda x: x.replace('-L', ""))
            merged_df = merged_df.pivot(index="filename", columns="acronym_Left", values=["Bilateral Density (cells/mm^3)"])
            merged_df.columns = merged_df.columns.get_level_values(1)
            merged_df.reset_index(drop=True, inplace=True)
            data_prewrangled[exp] = merged_df
        else:
            new_df = df[["filename","condition", "animal", "name", "acronym", "density (cells/mm^3)"]]
            new_df = new_df.pivot(index="filename", columns="acronym", values=["density (cells/mm^3)"])
            new_df.columns = new_df.columns.get_level_values(1)
            new_df.reset_index(drop=True, inplace=True)
            data_prewrangled[exp] = new_df
    if not ROI:
        data_wrangled = data_prewrangled
        print("Any ROIs.csv file have been found in the folder or you don't want anatomical organization, "
              "data won't be reorganized accordingly")
    else:
        data_wrangled = {exp:df[cols] for exp,df in data_prewrangled.items()}
    if ROI:
        if bilateral:
            pickle_save(data_wrangled, os.path.join(filepath, 'Data_Network_ROI_BL.pickle'))
        else:
            pickle_save(data_wrangled, os.path.join(filepath, 'Data_Network_ROI.pickle'))
    else:
        if bilateral:
            pickle_save(data_wrangled, os.path.join(filepath, 'Data_Network_BL.pickle'))
        else:
            pickle_save(data_wrangled,os.path.join(filepath, 'Data_Network.pickle'))
    for exp in data_wrangled.keys():
        if ROI:
            if bilateral:
                data_wrangled[exp].to_csv(os.path.join(filepath, str(exp) + '_Large_Network_ROI_BL.csv'), index = False)
            else:
                data_wrangled[exp].to_csv(os.path.join(filepath, str(exp) + '_Large_Network_ROI.csv'), index = False)

        else:
            if bilateral:
                data_wrangled[exp].to_csv(os.path.join(filepath, str(exp) + '_Large_Network_BL.csv'), index = False)
            else:
                data_wrangled[exp].to_csv(os.path.join(filepath, str(exp) + '_Large_Network.csv'), index = False)
    return data_wrangled

def dataManager(filepath = None, bilateral = False, ROI = False, processed = False):
    """
    Function that calls all the other functions for the data wrangling, also determines the
    filepath and output_folder

    Parameters
    ----------
    filepath: string containing the path were the csv files and the ROI file are located
    or None (the GUI would be used)
    Bilateral: Boolean that defines if you can to do bilateral analyses or all combine.
                True: separate left and right
                False: mean of left and right
    ROI: Boolean if you want to organize your data anatomically
    processed: Boolean that if you want to select the pickle of the data_wrangled for doing the network
                True: selection pickle
                False: data wrangling
    Returns
    -------
    dict or None
        A dictionary containing the data_wrangled dataframes per condition
        in your experiment (keys). All the data is saved in yoyr folder as a csv
        None is returned if any file has been selected.
    filepath: path of the files
    results_folder: path of a new generated folder with today's results

    """

    if not filepath:
        filepath = filefolderManager('folder', msg= ''.join(("Select the folder with the csv data from lifecanvas",
                                          "the discarded areas pickle (Discard.pickle') and the ROI.csv")))
        filepath = filepath[0] if not filepath == None and len(
            filepath) == 1 else sys.exit(''.join(('[Error] ',
                                                  'Only one output folder must be selected.')))
    results_folder = os.path.join(filepath, '_'.join((
        datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'),
        'network_analysis_plots')))
    os.mkdir(results_folder)
    if not processed:
        data_wrangled = data_wrangling(filepath, bilateral, ROI)
    else:
        data_wrangled = pickle_load(filefolderManager('file')[0])
    return data_wrangled, filepath, results_folder

"""
TEST
"""

filepath = input('Copy manually your filepath, if you want to use the GUI leave it blank')
data_wrangled, filepath, results_folder = dataManager(filepath, bilateral = False, ROI = False, processed =False)
