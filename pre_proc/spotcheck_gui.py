from tkinter import ttk
from tkinter import *
import tkinter
import os
from tkinter import filedialog


class GUI(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.grid()
        self.mainWindow = tkinter.Tk()
        self.mainWindow.title("Select LFP File")
        self.mainWindow.geometry('640x480-8-200')

        self.rightFrame = tkinter.Frame(self.mainWindow)
        self.rightFrame.grid(row=1, column=2, sticky='n')
        self.rightFrame.columnconfigure(0,weight=1)

        ttk.Button(self.rightFrame, text = "Select Directory", command=self.askdirectory).grid(row=0, column=0, padx=5, pady=5)

        # configure the columns
        self.mainWindow.columnconfigure(0, weight=1)
        self.mainWindow.columnconfigure(1, weight=1)
        self.mainWindow.grid_columnconfigure(2, weight=1)

        self.fileList = tkinter.Listbox(self.mainWindow)
        self.fileList.grid(row=1, column=0, sticky='nsew', rowspan=2)
        self.fileList.config(border=2, relief='sunken')


        self.listScroll = tkinter.Scrollbar(self.mainWindow, orient=tkinter.VERTICAL, command=self.fileList.yview)
        self.listScroll.grid(row=1, column=1, sticky='nsw', rowspan=2)
        self.fileList['yscrollcommand'] = self.listScroll.set

    # User selects directory with files to be analyzed
    def askdirectory(self):
        self.clear()
        self.directory = filedialog.askdirectory()
        print(self.directory)

        try:
            for zone in os.listdir(self.directory):
                self.fileList.insert(tkinter.END, zone)
        except FileNotFoundError:
            print('Try connecting to Server and try again')

    # Clear contents of the file list
    def clear(self):
        self.fileList.delete(0, 'end')

if __name__ == "__main__":
    guiFrame = GUI()
    guiFrame.mainloop()
