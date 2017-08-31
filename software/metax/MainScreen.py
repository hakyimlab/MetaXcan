__author__ = 'heroico'

import os
import logging
import weakref
import tkFileDialog
import tkMessageBox
import Tkinter
from subprocess import call
from threading import Thread
import metax.MainScreenView as MainScreenView
import metax.MetaXcanUITask as MetaXcanUITask
import metax.Exceptions as Exceptions
import MetaXcan
import metax.Formats as Formats
import metax.GWASUtilities as GWASUtilities
import metax.ZScoreCalculation as ZScoreCalculation
import metax.Normalization as Normalization
from Utilities import TS
from Utilities import checkSubdirectorySanity


GWAS_INPUT_DEFAULT = "data/GWAS"
BETA_FOLDER = "intermediate/beta"
COVARIANCE_FILE = "data/covariance.DGN-WB_0.5.txt.gz"
WEIGHT_DB_PATH = "data/DGN-WB_0.5.db"
OUTPUT_PATH = "results/zscores.csv"

SKIP = "Skip"

class MainScreen(object):
    def __init__(self, root, app):
        self.running = False
        self.root = root
        self.app = weakref.ref(app)

        self.process = None
        self.monitor = None
        self.poll = None

        self.cwd = os.getcwd()
        self.gwas_folder = "."
        if os.path.exists(GWAS_INPUT_DEFAULT):
            self.gwas_folder = GWAS_INPUT_DEFAULT

        self.beta_folder = BETA_FOLDER

        self.weight_db_path = "."
        if os.path.exists(WEIGHT_DB_PATH):
            self.weight_db_path = WEIGHT_DB_PATH

        self.output_path = OUTPUT_PATH

        # self.compressed_on = Tkinter.BooleanVar()
        # self.compressed_on.set(False)

        self.gwas_file_pattern_value = Tkinter.StringVar()

        self.separator_value = Tkinter.StringVar()

        self.snp_value = Tkinter.StringVar()
        self.snp_value.set("SNP")

        self.non_effect_allele_value = Tkinter.StringVar()
        self.non_effect_allele_value.set("A2")

        self.effect_allele_value = Tkinter.StringVar()
        self.effect_allele_value.set("A1")

        self.or_on = Tkinter.BooleanVar()
        self.or_on.set(False)
        self.or_value = Tkinter.StringVar()
        self.or_value.set("OR")

        self.beta_on = Tkinter.BooleanVar()
        self.beta_on.set(False)
        self.beta_value = Tkinter.StringVar()
        self.beta_value.set("BETA")

        self.beta_sign_on = Tkinter.BooleanVar()
        self.beta_sign_on.set(False)
        self.beta_sign_value = Tkinter.StringVar()
        self.beta_sign_value.set("BETA_SIGN")

        self.beta_z_on = Tkinter.BooleanVar()
        self.beta_z_on.set(False)
        self.beta_z_value = Tkinter.StringVar()
        self.beta_z_value.set("Z")

        self.p_on = Tkinter.BooleanVar()
        self.p_on.set(False)
        self.p_value = Tkinter.StringVar()
        self.p_value.set("P")

        self.se_on = Tkinter.BooleanVar()
        self.se_on.set(False)
        self.se_value = Tkinter.StringVar()
        self.se_value.set("SE")

        self.frequency_on = Tkinter.BooleanVar()
        self.frequency_on.set(True)
        self.frequency_value = Tkinter.StringVar()
        self.frequency_value.set("")

        self.covariance_file = "."
        if os.path.exists(COVARIANCE_FILE):
            self.covariance_file = COVARIANCE_FILE

        self.view = MainScreenView.MainScreenView(root, self)


###########
# UI events
###########
    def exclusiveToggle(self, boolean_var, entry, exclusive_vars, exclusive_entries):
        self.toggleOption(boolean_var, entry)
        on = boolean_var.get()
        for i,v in enumerate(exclusive_vars):
            e = exclusive_entries[i]
            if on:
                self.setOption(v, e, False)

    def setOption(self, boolean_var, entry, state):
        boolean_var.set(state)
        self.toggleOption(boolean_var, entry)

    def toggleOption(self, boolean_var, entry):
        on = boolean_var.get()
        entry.config(state=(Tkinter.NORMAL if on else Tkinter.DISABLED))

    def gwasFolderButtonPressed(self):
        rel = self.folderButtonPressed(self.view.gwas_button, True)
        if rel:
            self.gwas_folder = rel

    def betaFolderButtonPressed(self):
        rel = self.folderButtonPressed(self.view.beta_folder_button, False)
        if rel:
            self.beta_folder = rel

    def weightDBButtonPressed(self):
        rel = self.openFileButtonPressed(self.view.weight_db_button)
        if rel:
            self.weight_db_path = rel

    def outputButtonPressed(self):
        rel = self.saveFileButtonPressed(self.view.output_button)
        if rel:
            self.output_path = rel

    def covarianceFileButtonPressed(self):
        rel = self.openFileButtonPressed(self.view.covariance_file_button)
        if rel:
            self.covariance_file = rel

    def quitButtonPressed(self):
        self.root.quit()

    def actionButtonPressed(self):
        if not self.running:
            self.run()
        else:
            self.interrupt()

    def folderButtonPressed(self, button, must_exist):
        dir = tkFileDialog.askdirectory(mustexist=must_exist)
        if len(dir) == 0:
            return None
        rel = os.path.relpath(dir, self.cwd)
        button.config(text=rel)
        return rel

    def saveFileButtonPressed(self, button):
        file = tkFileDialog.asksaveasfilename()
        if len(file) == 0:
            return None
        rel = os.path.relpath(file, self.cwd)
        button.config(text=rel)
        return rel

    def openFileButtonPressed(self, button):
        file = tkFileDialog.askopenfilename()
        if len(file) == 0:
            return None
        rel = os.path.relpath(file, self.cwd)
        button.config(text=rel)
        return rel


#task
    def run(self):
        should_run = self.checkClearToRun()
        if not should_run:
            return
        self.running = True
        self.view.runMode()
        self.launchTask()

    def checkClearToRun(self):
        # sane = checkSubdirectorySanity(self.cwd, self.beta_folder)
        # if not sane:
        #     tkMessageBox.showwarning( "Beta Folder", "Beta folder cannot be current directory, or ancestor.")
        #     return False

        # clear, message, clean_up_beta, clean_up_results = self.checkGenerated()
        # if not clear:
        #     answer = tkMessageBox.askokcancel(TS("Warning!"), message, icon=tkMessageBox.ERROR)
        #     if answer:
        #         self.cleanUpGenerated(clean_up_beta, clean_up_results)
        #     else:
        #         return False

        return True

    def checkGenerated(self):
        clear = True

        beta_empty = True
        if os.path.exists(self.beta_folder):
            beta_empty = len(os.listdir(self.beta_folder)) == 0
        results_clear = not os.path.exists(self.output_path)

        clean_up_beta = False
        clean_up_results = False

        message = None
        if not beta_empty and not results_clear:
            clear = False
            clean_up_beta = True
            clean_up_results = True
            message = TS("Path for results already exists, and intermediate folder is already occupied."
                          "Should we delete them to move forward (potentially dangerous), or do you wish to cancel?")
        elif not results_clear:
            clear = False
            clean_up_results = True
            message = TS("Path for results already exists. Do you wish to cancel, or should we delete it to move forward (potentially dangerous)?")
        elif not beta_empty:
            clear = False
            clean_up_beta = True
            message = TS("Intermediate folder is already occupied. Should we delete it to move forward (potentially dangerous), or do you wish to cancel? (dangerous)")

        return clear, message, clean_up_beta, clean_up_results

    def cleanUpGenerated(self, clean_up_beta, clean_up_results):
        if clean_up_beta:
            command = "rm -rf " + self.beta_folder
            call(command.split())

        if clean_up_results:
            command = "rm -rf " + self.output_path
            call(command.split())

    def interrupt(self):
        logging.info("interrupting")
        if self.process and self.process.is_alive():
            self.process.terminate()
            while(True):
                code = self.process.exitcode
                if code != None and not self.process.is_alive():
                    break
        self.stop()

    def stop(self):
        if not self.running:
            logging.info("short stop")
            return

        logging.info("stopping")
        self.running = False
        self.view.configureMode()

        if self.monitor:
            self.monitor.give_up = True
            self.monitor = None

        self.process = None

    def launchTask(self):
        self.view.text.delete("1.0", Tkinter.END)
        self.running = True
        work = self.buildWork()
        process, monitor = MetaXcanUITask.runWork(work, self.taskCallback)
        self.process = process
        self.monitor = monitor

        poll = Thread(target=self.pollProcess)
        self.poll = poll
        poll.daemon = True
        poll.start()

    def taskCallback(self, queue):
        self.view.text.insert(Tkinter.END, queue.get())
        self.view.text.see(Tkinter.END)

    def pollProcess(self):
        code = None
        while(self.running):
            if self.process and self.process.is_alive():
                code = self.process.exitcode
                if code:
                    break
            else:
                break
        logging.info("poll stopping %s", str(code) if code else "-")
        self.stop()

# MetaxCan process
    def buildWork(self):
        class MetaXcanArgs(object):
            def __init__(self, source):
                self.verbosity = "10"
                self.weight_db_path = source.weight_db_path
                self.gwas_folder = source.gwas_folder

                self.snp_column = source.snp_value.get()
                self.non_effect_allele_column = source.non_effect_allele_value.get()
                self.effect_allele_column = source.effect_allele_value.get()

                self.or_column = source.or_value.get() if source.or_on.get() else None
                self.beta_column = source.beta_value.get() if source.beta_on.get() else None
                self.beta_sign_column = source.beta_sign_value.get() if source.beta_sign_on.get() else None
                self.zscore_column = source.beta_z_value.get() if source.beta_z_on.get() else None
                self.frequency_column = source.frequency_value.get() if source.frequency_on.get() else None
                self.se_column = source.se_value.get() if source.se_on.get() else None
                self.pvalue_column = source.p_value.get() if source.p_on.get() else None
                self.gwas_file_pattern = source.gwas_file_pattern_value.get() if len(source.gwas_file_pattern_value.get()) else None
                self.separator = source.separator_value.get() if len(source.separator_value.get()) else None

                # TODO: implement this
                self.skip_until_header = None
                self.throw = True

                self.verbosity = "10"
                self.remove_ens_version = False
                self.model_db_path = source.weight_db_path
                self.output_file = source.output_path
                self.covariance = source.covariance_file
                self.throw = True
                self.overwrite = True

        beta_args = MetaXcanArgs(source=self)

        #TODO: maybe connect stuff together so that M03 passes stuff to M04

        class MetaxcanWorkWrapper(object):
            def __init__(self,args):
                self.args = args

            def run(self):
                try:
                    MetaXcan.run(self.args)
                except Exceptions.ReportableException, e:
                    logging.error(e.msg)
                except Exception as e:
                    logging.info("Exception when running task: %s", str(e))
                finally:
                    pass

        class WorkWrapper(object):
            def __init__(self, works):
                self.works = works

            def run(self):
                try:
                    #delete as we go so that stuff gets deleted
                    self.works = list(reversed(self.works))
                    for i in xrange(len(self.works) - 1, -1, -1):
                        work = self.works[i]
                        work.run()
                        del self.works[i]
                except Exceptions.ReportableException, e:
                    logging.error(e.msg)
                except Exception as e:
                    logging.info("Exception when running task: %s", str(e))
                finally:
                    pass

        work = MetaxcanWorkWrapper(beta_args)
        return work