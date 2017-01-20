__author__ = 'heroico'

import weakref
from Tkinter import *
from ttk import *
from Utilities import TS

RHW = 70
LHW = 20

class MainScreenView(Frame):
    def __init__(self, root, controller):
        Frame.__init__(self, root)
        self.pack(fill=BOTH, expand=1)
        self.controller = weakref.ref(controller)
        self.root = root
        self.disabable = []
        self.disabable_status = {}
        self.initStyle()
        self.initUI()

    def addToDisabable(self, control):
        self.disabable.append(control)

    def initUI(self):
        self.initPreprocessUI()
        self.initSharedUI()
        self.initProcessUI()
        self.initBottomUI()

    def initStyle(self):
        self.style = Style()
        self.style.configure("Red.TButton", foreground="red")
        self.style.configure("RHW.TButton", background = "#fff", anchor=W, width=RHW)

    def initBottomUI(self):
        b_container = Frame(self, relief=RAISED, borderwidth=1)
        b_container.pack(fill=BOTH, expand=1)

        controller = self.controller()

        quit_button = Button(b_container, text=TS("QUIT"), command=controller.quitButtonPressed, style="Red.TButton")
        quit_button.pack(side=RIGHT)
        self.action_button = Button(b_container, text=TS("Run"), command=controller.actionButtonPressed)
        self.action_button.pack(side = LEFT)

        t_container = Frame(self, relief=RAISED, borderwidth=1)
        t_container.pack(fill=BOTH, expand=1)
        self.text = Text(t_container, height = 10)
        self.text.pack(fill=BOTH, side=TOP)

    def initPreprocessUI(self):
        self.preprocess_container = Frame(self, relief=RAISED, borderwidth=1)
        self.preprocess_container.pack(fill=BOTH, expand=1)

        options_frame = Frame(self.preprocess_container)
        options_frame.pack(fill=X, expand=1, side=TOP)

        self.initPreprocessOptionUI(options_frame)
        self.initPreprocessColumnUI(options_frame)
        self.initGWASFolderUI(options_frame)

    def initPreprocessOptionUI(self, frame):
        options_frame = Frame(frame)
        options_frame.pack(fill=X, expand=1, side=TOP)

        label = Label(options_frame, text=TS("GWAS File Options"), relief=RIDGE)
        label.pack(side=TOP,fill=X,expand=1)

        file_options_frame = Frame(options_frame)
        file_options_frame.pack(fill=X, expand=1, side=TOP)

        controller = self.controller()

        # compressed_on = controller.compressed_on
        # self.compressed_check = Checkbutton(file_options_frame, text=TS('Compressed GWAS'), variable=compressed_on)
        # self.compressed_check.grid(row=0, column=0, sticky=W+E)
        # self.addToDisabable(self.compressed_check)

        file_pattern_value = controller.gwas_file_pattern_value
        file_pattern_label = Label(file_options_frame, text=TS("GWAS Filename Pattern"))
        file_pattern_label.grid(row=1, column=0, sticky=W)
        self.file_pattern_entry = Entry(file_options_frame, textvariable=file_pattern_value)
        self.file_pattern_entry.grid(row=1, column=1, sticky=W+E)
        self.addToDisabable(self.file_pattern_entry)

        separator_value = controller.separator_value
        separator_label = Label(file_options_frame, text=TS("separator"))
        separator_label.grid(row=2, column=0, sticky=W)
        self.separator_entry = Entry(file_options_frame, textvariable=separator_value)
        self.separator_entry.grid(row=2, column=1, sticky=W+E)
        self.addToDisabable(self.separator_entry)

    def initPreprocessColumnUI(self, frame):
        preprocess_frame = Frame(frame, relief=RAISED, borderwidth=1)
        preprocess_frame.pack(fill=X, expand=1, side=TOP)

        label = Label(preprocess_frame, text=TS("Data Columns"), relief=RIDGE)
        label.pack(fill=X, expand=1, side=TOP)

        column_frame = Frame(preprocess_frame)
        column_frame.pack(fill=X, expand=1, side=TOP)

        controller = self.controller()

        #
        snp_label = Label(column_frame, text=TS("SNP id"))
        snp_label.grid(row=0, column=0, sticky=W)
        self.snp_entry = Entry(column_frame, textvariable=controller.snp_value)
        self.snp_entry.grid(row=0, column=1)
        self.addToDisabable(self.snp_entry)

        #
        non_effect_allele_label = Label(column_frame, text=TS("Non Effect Allele"))
        non_effect_allele_label.grid(row=1, column=0, sticky=W)
        self.non_effect_allele_entry = Entry(column_frame, textvariable=controller.non_effect_allele_value)
        self.non_effect_allele_entry.grid(row=1, column=1)
        self.addToDisabable(self.non_effect_allele_entry)

        #
        effect_allele_label = Label(column_frame, text=TS("Effect Allele"))
        effect_allele_label.grid(row=2, column=0, sticky=W)
        self.effect_allele_entry = Entry(column_frame, textvariable=controller.effect_allele_value)
        self.effect_allele_entry.grid(row=2, column=1)
        self.addToDisabable(self.effect_allele_entry)

        # the following three are mutually exclusive
        or_on = controller.or_on
        or_value = controller.or_value
        self.or_check = Checkbutton(column_frame, text=TS('Odd Ratios'), variable=or_on)
        self.or_check.grid(row=0, column=3, sticky=W+E)
        or_state = NORMAL if or_on.get() else DISABLED
        self.or_entry = Entry(column_frame, textvariable=or_value, state=or_state)
        self.or_entry.grid(row=0, column=4, sticky=W)
        self.addToDisabable(self.or_check)
        self.addToDisabable(self.or_entry)

        beta_on = controller.beta_on
        beta_value = controller.beta_value
        self.beta_check = Checkbutton(column_frame, text=TS('Beta'), variable=beta_on)
        self.beta_check.grid(row=1, column=3, sticky=W+E)
        beta_state = NORMAL if beta_on.get() else DISABLED
        self.beta_entry = Entry(column_frame, textvariable=beta_value, state=beta_state)
        self.beta_entry.grid(row=1, column=4, sticky=W)
        self.addToDisabable(self.beta_check)
        self.addToDisabable(self.beta_entry)

        beta_sign_on = controller.beta_sign_on
        beta_sign_value = controller.beta_sign_value
        self.beta_sign_check = Checkbutton(column_frame, text=TS('Beta Sign'), variable=beta_sign_on)
        self.beta_sign_check.grid(row=2, column=3, sticky=W+E)
        beta_sign_state = NORMAL if beta_sign_on.get() else DISABLED
        self.beta_sign_entry = Entry(column_frame, textvariable=beta_sign_value, state=beta_sign_state)
        self.beta_sign_entry.grid(row=2, column=4, sticky=W)
        self.addToDisabable(self.beta_sign_check)
        self.addToDisabable(self.beta_sign_entry)

        or_cmd = lambda: controller.exclusiveToggle(or_on, self.or_entry, [beta_on, beta_sign_on], [self.beta_entry, self.beta_sign_entry])
        self.or_check.config(command=or_cmd)

        beta_lambda = lambda: controller.exclusiveToggle(beta_on, self.beta_entry, [or_on, beta_sign_on], [self.or_entry, self.beta_sign_entry])
        self.beta_check.config(command=beta_lambda)

        beta_sign_lambda = lambda: controller.exclusiveToggle(beta_sign_on, self.beta_sign_entry, [beta_on, or_on], [self.or_entry, self.beta_entry])
        self.beta_sign_check.config(command=beta_sign_lambda)

        #
        self.p_check = Checkbutton(column_frame, text=TS('P Value'), variable=controller.p_on)
        self.p_check.grid(row=0, column=5, sticky=W+E)
        p_state = NORMAL if controller.p_on.get() else DISABLED
        self.p_entry = Entry(column_frame, textvariable=controller.p_value, state = p_state)
        self.p_entry.grid(row=0, column=6, sticky=W)
        p_lambda = lambda : controller.toggleOption(controller.p_on, self.p_entry)
        self.p_check.config(command=p_lambda)
        self.addToDisabable(self.p_check)
        self.addToDisabable(self.p_entry)

    def initGWASFolderUI(self, options_frame):
        gwas_frame = Frame(options_frame)
        gwas_frame.pack(fill=X, expand=1, side=TOP)

        gwas_label = Label(gwas_frame, text=TS("GWAS data folder"), anchor=W)
        gwas_label.pack(side=LEFT)

        controller = self.controller()
        self.gwas_button = Button(gwas_frame, text=controller.gwas_folder,
            command=controller.gwasFolderButtonPressed,
            style="RHW.TButton"
        )
        self.gwas_button.pack(side=RIGHT, fill=X)
        self.addToDisabable(self.gwas_button)

    def initSharedUI(self):
        shared_container = Frame(self, relief=RAISED, borderwidth=1)
        shared_container.pack(fill=BOTH, expand=1)

        #self.initBetaFolderUI(shared_container)
        self.initWeightDBUI(shared_container)

    def initBetaFolderUI(self, options_frame):
        output_frame = Frame(options_frame)
        output_frame.pack(fill=X, expand=1, side=TOP)

        controller = self.controller()

        beta_label = Label(output_frame, text=TS("Processed beta folder"))
        beta_label.pack(side=LEFT)

        self.beta_folder_button = Button(output_frame, text=controller.beta_folder,
            command=controller.betaFolderButtonPressed,
            style="RHW.TButton"
        )
        self.beta_folder_button.pack(side=RIGHT)
        self.addToDisabable(self.beta_folder_button)

    def initWeightDBUI(self, frame):
        w_frame = Frame(frame)
        w_frame.pack(fill=X, expand=1, side=TOP)

        controller = self.controller()

        label = Label(w_frame, text=TS("Weight Model DB"))
        label.pack(side=LEFT)

        self.weight_db_button = Button(w_frame, text=controller.weight_db_path,
            command=controller.weightDBButtonPressed,
            style="RHW.TButton"
        )
        self.weight_db_button.pack(side=RIGHT)
        self.addToDisabable(self.weight_db_button)

    def initProcessUI(self):
        process_container = Frame(self, relief=RAISED, borderwidth=1)
        process_container.pack(fill=BOTH, expand=1)

        self.initCovarianceUI(process_container)
        self.initOutputFileUI(process_container)

    def initOutputFileUI(self, frame):
        output_frame = Frame(frame)
        output_frame.pack(fill=X, expand=1, side=TOP)
        controller = self.controller()

        label = Label(output_frame, text=TS("Output file"))
        label.pack(side=LEFT)

        self.output_button = Button(output_frame, text=controller.output_path,
            command = controller.outputButtonPressed,
            style = "RHW.TButton",
        )
        self.output_button.pack(side=RIGHT)
        self.addToDisabable(self.output_button)

    def initCovarianceUI(self, frame):
        output_frame = Frame(frame)
        output_frame.pack(fill=X, expand=1, side=TOP)
        controller = self.controller()

        label = Label(output_frame, text=TS("Covariance File"))
        label.pack(side=LEFT)

        self.covariance_file_button = Button(output_frame, text=controller.covariance_file,
            command = controller.covarianceFileButtonPressed,
            style = "RHW.TButton",
        )
        self.covariance_file_button.pack(side=RIGHT)
        self.addToDisabable(self.covariance_file_button)

    def configureMode(self):
        text = TS("Run")
        self.action_button.config(text=text)
        for control in self.disabable:
            state = self.disabable_status[str(control)]
            control.config(state=state)

    def runMode(self):
        text = TS("Stop")
        self.action_button.config(text=text)
        for control in self.disabable:
            self.disabable_status[str(control)] = control.cget("state")
            control.config(state=DISABLED)

