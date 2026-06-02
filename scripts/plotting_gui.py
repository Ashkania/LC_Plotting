#!/usr/bin/env python3
import os
import sys
import shlex
import threading
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext


def resource_path(relpath: str) -> str:
    # Return path relative to this script (scripts/)
    return os.path.join(os.path.dirname(__file__), relpath)


class PlottingGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('LC Plotting GUI')
        self.geometry('820x680')
        self.create_widgets()

    def create_widgets(self):
        frm = ttk.Frame(self)
        frm.pack(fill='both', expand=True, padx=8, pady=8)

        # File selection
        file_row = ttk.Frame(frm)
        file_row.pack(fill='x', pady=4)
        ttk.Label(file_row, text='HDF5 file:').pack(side='left')
        self.file_var = tk.StringVar()
        ttk.Entry(file_row, textvariable=self.file_var, width=70).pack(side='left', padx=6)
        ttk.Button(file_row, text='Browse', command=self.browse_file).pack(side='left')

        opts_row = ttk.Frame(frm)
        opts_row.pack(fill='x', pady=4)

        # Aperture
        ttk.Label(opts_row, text='Aperture:').grid(row=0, column=0, sticky='w')
        self.aperture_var = tk.IntVar(value=4)
        ttk.Spinbox(opts_row, from_=0, to=20, textvariable=self.aperture_var, width=6).grid(row=0, column=1)

        # Period & epoch
        ttk.Label(opts_row, text='Period:').grid(row=0, column=2, sticky='w')
        self.period_var = tk.StringVar()
        ttk.Entry(opts_row, textvariable=self.period_var, width=12).grid(row=0, column=3)
        ttk.Label(opts_row, text='Epoch:').grid(row=0, column=4, sticky='w')
        self.epoch_var = tk.StringVar()
        ttk.Entry(opts_row, textvariable=self.epoch_var, width=12).grid(row=0, column=5)

        # Mag types checkboxes
        mags_row = ttk.LabelFrame(frm, text='Magnitude types')
        mags_row.pack(fill='x', pady=6)
        self.magfit_var = tk.BooleanVar(value=True)
        self.epd_var = tk.BooleanVar(value=False)
        self.tfa_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(mags_row, text='magfit', variable=self.magfit_var).pack(side='left', padx=6)
        ttk.Checkbutton(mags_row, text='epd', variable=self.epd_var).pack(side='left', padx=6)
        ttk.Checkbutton(mags_row, text='tfa', variable=self.tfa_var).pack(side='left', padx=6)

        # sep-by and mode
        sel_row = ttk.Frame(frm)
        sel_row.pack(fill='x', pady=4)
        ttk.Label(sel_row, text='Separate by:').pack(side='left')
        self.sep_var = tk.StringVar(value='photref')
        ttk.Radiobutton(sel_row, text='photref', variable=self.sep_var, value='photref').pack(side='left')
        ttk.Radiobutton(sel_row, text='channel', variable=self.sep_var, value='channel').pack(side='left')

        ttk.Label(sel_row, text='Mode:').pack(side='left', padx=(20, 0))
        self.mode_var = tk.StringVar(value='single')
        ttk.Radiobutton(sel_row, text='single', variable=self.mode_var, value='single').pack(side='left')
        ttk.Radiobutton(sel_row, text='folded', variable=self.mode_var, value='folded').pack(side='left')

        # Binning
        bin_row = ttk.Frame(frm)
        bin_row.pack(fill='x', pady=4)
        ttk.Label(bin_row, text='Binning:').pack(side='left')
        self.binning_method = tk.StringVar(value='None')
        ttk.Combobox(bin_row, textvariable=self.binning_method, values=['None', 'time', 'phase'], width=8).pack(side='left', padx=6)
        ttk.Label(bin_row, text='Size:').pack(side='left')
        self.binning_size = tk.StringVar()
        ttk.Entry(bin_row, textvariable=self.binning_size, width=8).pack(side='left')

        # selected indices
        selind_row = ttk.Frame(frm)
        selind_row.pack(fill='x', pady=4)
        ttk.Label(selind_row, text='Selected photref/channel indices (space/comma separated):').pack(side='left')
        self.selected_var = tk.StringVar()
        ttk.Entry(selind_row, textvariable=self.selected_var, width=30).pack(side='left', padx=6)

        # Flags
        flags_row = ttk.Frame(frm)
        flags_row.pack(fill='x', pady=6)
        self.tess_var = tk.BooleanVar()
        self.gaia_var = tk.BooleanVar()
        self.save_var = tk.StringVar(value='show')
        ttk.Checkbutton(flags_row, text='TESS plot', variable=self.tess_var).pack(side='left', padx=6)
        ttk.Checkbutton(flags_row, text='GAIA plot', variable=self.gaia_var).pack(side='left', padx=6)
        ttk.Radiobutton(flags_row, text='Show', variable=self.save_var, value='show').pack(side='left', padx=6)
        ttk.Radiobutton(flags_row, text='Save', variable=self.save_var, value='save').pack(side='left', padx=6)

        # Y limits
        ylim_row = ttk.Frame(frm)
        ylim_row.pack(fill='x', pady=4)
        ttk.Label(ylim_row, text='Ylim (Ymin Ymax):').pack(side='left')
        self.ylim_min = tk.StringVar()
        self.ylim_max = tk.StringVar()
        ttk.Entry(ylim_row, textvariable=self.ylim_min, width=10).pack(side='left', padx=4)
        ttk.Entry(ylim_row, textvariable=self.ylim_max, width=10).pack(side='left', padx=4)

        # Periodic text
        text_row = ttk.Frame(frm)
        text_row.pack(fill='x', pady=4)
        ttk.Label(text_row, text='Text (optional):').pack(side='left')
        self.text_var = tk.StringVar()
        ttk.Entry(text_row, textvariable=self.text_var, width=50).pack(side='left', padx=6)

        # Run and output
        run_row = ttk.Frame(frm)
        run_row.pack(fill='x', pady=8)
        ttk.Button(run_row, text='Run', command=self.run_plot).pack(side='left')
        ttk.Button(run_row, text='Clear output', command=self.clear_output).pack(side='left', padx=6)

        self.output = scrolledtext.ScrolledText(frm, height=20)
        self.output.pack(fill='both', expand=True, pady=6)

    def browse_file(self):
        path = filedialog.askopenfilename(title='Select HDF5 file', filetypes=[('HDF5', '*.h5 *.hdf5 *.fits.fz'), ('All files', '*.*')])
        if path:
            self.file_var.set(path)

    def clear_output(self):
        self.output.delete('1.0', tk.END)

    def append_output(self, text):
        self.output.insert(tk.END, text)
        self.output.see(tk.END)

    def build_args(self):
        path = self.file_var.get().strip()
        if not path:
            messagebox.showerror('Missing file', 'Please select an HDF5 file')
            return None

        args = [sys.executable, resource_path('plotting.py'), path]
        # mag types
        mag_types = []
        if self.magfit_var.get():
            mag_types.append('magfit')
        if self.epd_var.get():
            mag_types.append('epd')
        if self.tfa_var.get():
            mag_types.append('tfa')
        if mag_types:
            args.append('--mag-types')
            args.extend(mag_types)

        args += ['--aperture', str(self.aperture_var.get())]

        if self.period_var.get().strip():
            args += ['--period', self.period_var.get().strip()]
        if self.epoch_var.get().strip():
            args += ['--epoch', self.epoch_var.get().strip()]

        if self.text_var.get().strip():
            args += ['--text', self.text_var.get().strip()]

        args += ['--sep-by', self.sep_var.get()]
        args += ['--mode', self.mode_var.get()]

        if self.binning_method.get() != 'None' and self.binning_size.get().strip():
            args += ['--binning', self.binning_method.get(), self.binning_size.get().strip()]

        if self.selected_var.get().strip():
            # split by comma/space
            tokens = [t for t in shlex.split(self.selected_var.get().replace(',', ' '))]
            if tokens:
                args.append('--selected')
                args.extend(tokens)

        if self.tess_var.get():
            args.append('--TESS-plot')
        if self.gaia_var.get():
            args.append('--GAIA-plot')

        args += ['--save-or-show-plot', self.save_var.get()]

        if self.ylim_min.get().strip() and self.ylim_max.get().strip():
            args += ['--ylim', self.ylim_min.get().strip(), self.ylim_max.get().strip()]

        return args

    def run_plot(self):
        args = self.build_args()
        if args is None:
            return

        # Run in a thread to keep UI responsive
        def target():
            self.append_output('Running: ' + ' '.join(shlex.quote(a) for a in args) + '\n')
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in proc.stdout:
                self.append_output(line)
            proc.wait()
            self.append_output(f"\nProcess exited with code {proc.returncode}\n")

        threading.Thread(target=target, daemon=True).start()


def main():
    app = PlottingGUI()
    app.mainloop()


if __name__ == '__main__':
    main()
