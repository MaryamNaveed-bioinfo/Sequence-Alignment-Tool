import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
from msa_tool import SequenceAlignmentTool

# ── Colours ──────────────────────────────────────────────────────
BG = "#0d1117"
PANEL = "#161b22"
CARD = "#1c2330"
INP = "#0d1520"
WIRE = "#21303f"
GLOW = "#00ff9d"
CYAN = "#00e5ff"
AMBER = "#ffb700"
RED = "#ff4d6d"
TXT = "#c8f0e0"
DIM = "#4a7060"

NUC = {"A": "#00ff9d", "T": "#ff4d6d", "U": "#ff4d6d",
       "G": "#ffb700", "C": "#00e5ff", "N": "#556655"}

FM = ("Courier New", 11)
FMB = ("Courier New", 11, "bold")
FMS = ("Courier New", 10)
FUI = ("Courier New",  9)

SCORE_TIPS = {
    "match":      "+reward when bases match.  DNA: 1-2  Protein: 1-3",
    "mismatch":   "Penalty for different bases.  DNA: -1 to -2",
    "gap_open":   "Penalty to OPEN a new gap (paid once).  EMBOSS default: -10",
    "gap_extend": "Penalty per position EXTENDED in gap.  EMBOSS default: -0.5",
}


# ── Tooltip ───────────────────────────────────────────────────────
class Tip(tk.Toplevel):
    _i = None

    @classmethod
    def show(cls, w, txt):
        cls.hide()
        t = tk.Toplevel(w)
        t.wm_overrideredirect(True)
        t.wm_geometry(
            f"+{w.winfo_rootx()+w.winfo_width()+6}+{w.winfo_rooty()}")
        tk.Label(t, text=txt, font=FMS, bg=CARD, fg=CYAN,
                 padx=10, pady=6, justify="left").pack(padx=1, pady=1)
        cls._i = t

    @classmethod
    def hide(cls):
        if cls._i:
            try:
                cls._i.destroy()
            except Exception:
                pass
            cls._i = None


# ── Main App ──────────────────────────────────────────────────────
class MSAApp:
    def __init__(self, root):
        self.root = root
        root.title("Sequence Alignment Tool")
        root.geometry("1480x940")
        root.configure(bg=BG)
        root.minsize(1100, 720)

        self.aligner = SequenceAlignmentTool()
        self.seq_boxes = []
        self.last_out = ""

        self._style()
        self._menu()
        self._build()

    # ── TTK style ────────────────────────────────────────────────
    def _style(self):
        s = ttk.Style()
        s.theme_use("clam")
        s.configure("TRadiobutton", background=PANEL, foreground=TXT,
                    font=FUI, focuscolor=PANEL)
        s.map("TRadiobutton", background=[("active", PANEL)],
              foreground=[("active", GLOW)])
        s.configure("TNotebook", background=BG, borderwidth=0)
        s.configure("TNotebook.Tab", background=CARD, foreground=DIM,
                    font=FUI, padding=(12, 6))
        s.map("TNotebook.Tab", background=[("selected", BG)],
              foreground=[("selected", GLOW)])
        s.configure("Treeview", background=CARD, foreground=TXT,
                    fieldbackground=CARD, rowheight=22, font=FMS)
        s.configure("Treeview.Heading", background=PANEL, foreground=CYAN,
                    font=("Courier New", 9, "bold"), relief="flat")
        s.map("Treeview", background=[("selected", WIRE)])
        s.configure("TProgressbar", troughcolor=CARD, background=GLOW)
        s.configure("TScrollbar", background=PANEL, troughcolor=BG,
                    arrowcolor=DIM, relief="flat")

    # ── Menu ─────────────────────────────────────────────────────
    def _menu(self):
        mb = tk.Menu(self.root, bg=PANEL, fg=TXT,
                     activebackground=WIRE, activeforeground=GLOW, relief="flat")
        self.root.config(menu=mb)
        fm = tk.Menu(mb, tearoff=0, bg=PANEL, fg=TXT,
                     activebackground=WIRE, activeforeground=GLOW)
        fm.add_command(label="Open FASTA",
                       command=self._load_fasta, accelerator="Ctrl+O")
        fm.add_command(label="Save Output", command=self._save,
                       accelerator="Ctrl+S")
        fm.add_separator()
        fm.add_command(label="Exit", command=self.root.quit)
        mb.add_cascade(label="File", menu=fm)
        self.root.bind("<Control-o>", lambda e: self._load_fasta())
        self.root.bind("<Control-s>", lambda e: self._save())

    # ── Layout ───────────────────────────────────────────────────
    def _build(self):
        # Header
        hdr = tk.Frame(self.root, bg=PANEL, height=46)
        hdr.pack(fill="x")
        hdr.pack_propagate(False)
        tk.Frame(hdr, bg=GLOW, width=3).pack(side="left", fill="y")
        tk.Label(hdr, text="  SEQUENCE ALIGNMENT TOOL",
                 font=("Courier New", 14, "bold"), bg=PANEL, fg=GLOW).pack(side="left", padx=8)
        tk.Label(hdr, text="NW · SW · MSA  |  EMBOSS-Compatible  |  MEGA Colors",
                 font=FUI, bg=PANEL, fg=DIM).pack(side="left")
        cf = tk.Frame(hdr, bg=PANEL)
        cf.pack(side="right", padx=14)
        for b, c in [("A", NUC["A"]), ("T", NUC["T"]), ("G", NUC["G"]), ("C", NUC["C"])]:
            tk.Label(cf, text=f" {b} ", font=FMB, bg=c, fg=BG).pack(
                side="left", padx=2, pady=10)
        tk.Frame(self.root, bg=GLOW, height=1).pack(fill="x")

        # Body
        body = tk.PanedWindow(self.root, orient="horizontal", bg=BG,
                              sashwidth=4, sashrelief="flat")
        body.pack(fill="both", expand=True)

        left = tk.Frame(body, bg=PANEL, width=440)
        body.add(left, minsize=400)
        self._left_panel(left)

        right = tk.Frame(body, bg=BG)
        body.add(right, minsize=600)
        self._output_panel(right)

        # Status bar
        tk.Frame(self.root, bg=WIRE, height=1).pack(fill="x")
        sb = tk.Frame(self.root, bg=PANEL, height=24)
        sb.pack(fill="x", side="bottom")
        sb.pack_propagate(False)
        self._status_lbl = tk.Label(sb, text="Ready", font=FUI,
                                    bg=PANEL, fg=DIM, anchor="w")
        self._status_lbl.pack(side="left", padx=10, fill="y")
        self._pb = ttk.Progressbar(sb, mode="indeterminate",
                                   length=100, style="TProgressbar")

    # ── Left panel (scrollable) ───────────────────────────────────
    def _left_panel(self, parent):
        c = tk.Canvas(parent, bg=PANEL, bd=0, highlightthickness=0)
        sb = ttk.Scrollbar(parent, orient="vertical", command=c.yview)
        c.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        c.pack(fill="both", expand=True)
        inner = tk.Frame(c, bg=PANEL)
        cw = c.create_window((0, 0), window=inner, anchor="nw")
        c.bind("<Configure>", lambda e: c.itemconfig(cw, width=e.width))
        inner.bind("<Configure>", lambda e: c.configure(
            scrollregion=c.bbox("all")))
        c.bind_all("<MouseWheel>", lambda e: c.yview_scroll(
            int(-e.delta/120), "units"))

        self._section_algo(inner)
        self._section_params(inner)
        self._section_input(inner)
        self._section_buttons(inner)
        self._section_stats(inner)

    def _sec(self, parent, title, color=GLOW):
        tk.Label(parent, text=title, font=("Courier New", 8, "bold"),
                 bg=PANEL, fg=color).pack(anchor="w", padx=10, pady=(12, 0))
        tk.Frame(parent, bg=color, height=1).pack(
            fill="x", padx=10, pady=(1, 4))
        card = tk.Frame(parent, bg=CARD, padx=8, pady=6,
                        highlightthickness=1, highlightbackground=WIRE)
        card.pack(fill="x", padx=8, pady=(0, 4))
        return card

    # ── Algorithm ────────────────────────────────────────────────
    def _section_algo(self, p):
        body = self._sec(p, "ALGORITHM", GLOW)
        self.aln_mode = tk.StringVar(value="global")
        for val, lbl in [("global", "Global — Needleman-Wunsch"),
                         ("local",  "Local  — Smith-Waterman"),
                         ("msa",    "Multi  — Progressive MSA")]:
            ttk.Radiobutton(body, text=lbl, variable=self.aln_mode,
                            value=val, command=self._rebuild_input).pack(anchor="w", pady=2)

    # ── Scoring params ───────────────────────────────────────────
    def _section_params(self, p):
        body = self._sec(p, "SCORING PARAMETERS", CYAN)
        self._pvars = {}
        params = [
            ("match",      "Match Score",      2),
            ("mismatch",   "Mismatch Penalty", -1),
            ("gap_open",   "Gap Open",         -10),
            ("gap_extend", "Gap Extend",       -0.5),
        ]
        for key, lbl, default in params:
            row = tk.Frame(body, bg=CARD)
            row.pack(fill="x", pady=3)
            tk.Label(row, text=lbl, font=FUI, bg=CARD, fg=TXT,
                     width=16, anchor="w").pack(side="left")
            var = tk.StringVar(value=str(default))
            e = tk.Entry(row, textvariable=var, width=6, font=FMB,
                         bg=INP, fg=AMBER, insertbackground=GLOW,
                         relief="flat", bd=3, justify="center",
                         highlightthickness=1, highlightbackground=WIRE,
                         highlightcolor=GLOW)
            e.pack(side="left", padx=6)
            tip = tk.Label(row, text="?", font=FUI,
                           bg=CARD, fg=CYAN, cursor="hand2")
            tip.pack(side="left")
            tip.bind("<Enter>", lambda e, k=key,
                     w=tip: Tip.show(w, SCORE_TIPS[k]))
            tip.bind("<Leave>", lambda e: Tip.hide())
            self._pvars[key] = var

        # Presets
        pr = tk.Frame(body, bg=CARD)
        pr.pack(fill="x", pady=(6, 0))
        tk.Label(pr, text="Presets:", font=FUI, bg=CARD,
                 fg=DIM).pack(side="left", padx=(0, 6))
        # (label, match, mismatch, gap_open, gap_extend)
        presets = [
            ("EMBOSS DNA",     2,  -1,   -10, -0.5),
            ("EMBOSS Strict",  1,  -3,   -10, -0.5),
            ("Protein",        3,  -1,   -12, -1.0),
            ("Fast/Simple",    2,  -1,    -2, -0.5),
        ]
        for lbl, m, mm, go, ge in presets:
            b = tk.Label(pr, text=f"[{lbl}]", font=FUI,
                         bg=CARD, fg=CYAN, cursor="hand2")
            b.pack(side="left", padx=3)
            b.bind("<Button-1>", lambda e, m=m, mm=mm,
                   go=go, ge=ge: self._preset(m, mm, go, ge))
            b.bind("<Enter>", lambda e, b=b: b.config(fg=GLOW))
            b.bind("<Leave>", lambda e, b=b: b.config(fg=CYAN))

    def _preset(self, m, mm, go, ge):
        self._pvars["match"].set(str(m))
        self._pvars["mismatch"].set(str(mm))
        self._pvars["gap_open"].set(str(go))
        self._pvars["gap_extend"].set(str(ge))

    # ── Sequence input ───────────────────────────────────────────
    def _section_input(self, p):
        outer = self._sec(p, "SEQUENCE INPUT", AMBER)

        # --- Static part: file load button (never destroyed) ---
        frow = tk.Frame(outer, bg=CARD)
        frow.pack(fill="x", pady=(0, 4))
        tk.Button(frow, text="Load FASTA File", command=self._load_fasta,
                  bg=AMBER, fg=BG, font=FMB, relief="flat", bd=0,
                  padx=12, pady=5, cursor="hand2",
                  activebackground=GLOW, activeforeground=BG).pack(side="left")
        tk.Label(frow, text="  or paste below",
                 font=FUI, bg=CARD, fg=DIM).pack(side="left")
        tk.Frame(outer, bg=WIRE, height=1).pack(fill="x", pady=(0, 6))

        # --- Dynamic part: text boxes (rebuilt on mode change) ---
        self._dyn_area = tk.Frame(outer, bg=CARD)
        self._dyn_area.pack(fill="x")

        self.seq_boxes = []
        self._rebuild_input()

    def _rebuild_input(self):
        # Only destroy/recreate the DYNAMIC text-box area, not the file button
        for w in self._dyn_area.winfo_children():
            w.destroy()
        self.seq_boxes = []
        mode = self.aln_mode.get()
        if mode in ("global", "local"):
            self._mkbox(self._dyn_area, "Sequence 1", 5)
            self._mkbox(self._dyn_area, "Sequence 2", 5)
        else:
            self._mkbox(self._dyn_area,
                        "FASTA input  (≥ 3 sequences, FASTA format)", 14)

    def _mkbox(self, parent, label, rows):
        tk.Label(parent, text=label, font=FUI,
                 bg=CARD, fg=DIM).pack(anchor="w")
        bdr = tk.Frame(parent, bg=WIRE, padx=1, pady=1)
        bdr.pack(fill="x", pady=(1, 8))
        xsb = ttk.Scrollbar(bdr, orient="horizontal")
        box = tk.Text(bdr, height=rows, font=FM, bg=INP, fg=GLOW,
                      insertbackground=GLOW, selectbackground=WIRE,
                      relief="flat", bd=3, wrap="none",
                      highlightthickness=0,
                      xscrollcommand=xsb.set)
        xsb.configure(command=box.xview)
        box.pack(fill="x")
        xsb.pack(fill="x")
        box.bind("<FocusIn>", lambda e, f=bdr: f.config(bg=GLOW))
        box.bind("<FocusOut>", lambda e, f=bdr: f.config(bg=WIRE))
        self.seq_boxes.append(box)

    # ── Buttons ──────────────────────────────────────────────────
    def _section_buttons(self, p):
        body = self._sec(p, "ACTIONS", GLOW)
        row = tk.Frame(body, bg=CARD)
        row.pack(fill="x")
        for txt, cmd, col in [("▶ Run",     self._run_threaded, GLOW),
                              ("⟳ Clear",  self._clear,        CYAN),
                              ("↓ Export", self._save,         AMBER)]:
            b = tk.Button(row, text=txt, command=cmd,
                          bg=col, fg=BG, font=FMB,
                          relief="flat", bd=0, padx=14, pady=6,
                          cursor="hand2", activebackground=col, activeforeground=BG)
            b.pack(side="left", padx=(0, 6))
            b.bind("<Enter>", lambda e, b=b,
                   c=col: b.config(bg=self._brighten(c)))
            b.bind("<Leave>", lambda e, b=b, c=col: b.config(bg=c))

        opt = tk.Frame(body, bg=CARD)
        opt.pack(fill="x", pady=(8, 0))
        tk.Label(opt, text="Color:", font=FUI,
                 bg=CARD, fg=DIM).pack(side="left")
        self.cmode = tk.StringVar(value="nucleotide")
        for val, lbl in [("nucleotide", "MEGA NUC"), ("matchmiss", "Match/Mismatch")]:
            ttk.Radiobutton(opt, text=lbl, variable=self.cmode,
                            value=val).pack(side="left", padx=4)
        tk.Label(opt, text="  Block:", font=FUI,
                 bg=CARD, fg=DIM).pack(side="left")
        self.chunk = tk.IntVar(value=60)
        tk.Spinbox(opt, from_=20, to=120, textvariable=self.chunk,
                   width=4, font=FUI, bg=INP, fg=AMBER,
                   buttonbackground=CARD, relief="flat",
                   highlightthickness=0).pack(side="left", padx=4)

    def _brighten(self, hex_col):
        r, g, b = int(hex_col[1:3], 16), int(
            hex_col[3:5], 16), int(hex_col[5:7], 16)
        r, g, b = min(255, r+30), min(255, g+30), min(255, b+30)
        return f"#{r:02x}{g:02x}{b:02x}"

    # ── Stats ────────────────────────────────────────────────────
    def _section_stats(self, p):
        body = self._sec(p, "STATISTICS  (EMBOSS-style)", CYAN)
        self._slbls = {}
        grid = tk.Frame(body, bg=CARD)
        grid.pack(fill="x")
        items = [
            ("Score",      "score"),
            ("Type",       "seq_type"),
            ("Aln Length", "length"),
            ("Identity",   "identity"),
            ("Similarity", "similarity"),
            ("Matches",    "matches"),
            ("Mismatches", "mismatches"),
            ("Gaps",       "gaps"),
            ("Gap %",      "gap_pct"),
        ]
        for i, (lbl, key) in enumerate(items):
            r, c = divmod(i, 2)
            c0 = c * 2
            tk.Label(grid, text=lbl+":", font=FUI, bg=CARD, fg=DIM,
                     anchor="e", width=11).grid(row=r, column=c0, padx=(4, 2), pady=2, sticky="e")
            vl = tk.Label(grid, text="—", font=FMB, bg=CARD,
                          fg=GLOW, anchor="w", width=11)
            vl.grid(row=r, column=c0+1, padx=(0, 6), pady=2, sticky="w")
            self._slbls[key] = vl

        # EMBOSS note
        tk.Label(body, text="Identity = matches / aln_length  (matches EMBOSS)",
                 font=("Courier New", 7), bg=CARD, fg=DIM).pack(anchor="w", padx=4, pady=(4, 0))

    def _upd_stats(self, d):
        self._slbls["score"].config(text=str(d.get("score", "—")))
        self._slbls["seq_type"].config(text=d.get("seq_type", "—"))
        self._slbls["length"].config(text=str(d.get("length", "—")))
        self._slbls["identity"].config(text=f"{d.get('identity', 0):.2f}%")
        self._slbls["similarity"].config(text=f"{d.get('similarity', 0):.2f}%")
        self._slbls["matches"].config(text=str(d.get("matches", "—")))
        self._slbls["mismatches"].config(text=str(d.get("mismatches", "—")))
        self._slbls["gaps"].config(text=str(d.get("gaps", "—")))
        self._slbls["gap_pct"].config(text=f"{d.get('gap_percent', 0):.1f}%")

    # ── Output panel ─────────────────────────────────────────────
    def _output_panel(self, parent):
        nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True, padx=4, pady=4)

        t1 = tk.Frame(nb, bg=BG)
        nb.add(t1, text="  Alignment  ")
        t2 = tk.Frame(nb, bg=BG)
        nb.add(t2, text="  Score Matrix  ")
        t3 = tk.Frame(nb, bg=BG)
        nb.add(t3, text="  Seq Info  ")

        # Alignment tab
        self.out = tk.Text(t1, wrap="none", font=FM, bg=BG, fg=TXT,
                           insertbackground=GLOW, selectbackground=WIRE,
                           relief="flat", bd=0, highlightthickness=0, state="disabled")
        xsb = ttk.Scrollbar(t1, orient="horizontal", command=self.out.xview)
        ysb = ttk.Scrollbar(t1, orient="vertical",   command=self.out.yview)
        self.out.configure(xscrollcommand=xsb.set, yscrollcommand=ysb.set)
        xsb.pack(side="bottom", fill="x")
        ysb.pack(side="right",  fill="y")
        self.out.pack(fill="both", expand=True, padx=4, pady=4)
        self._setup_tags()

        # Matrix tab
        cols = ("Seq 1", "Seq 2", "Score", "Identity",
                "Similarity", "Matches", "Mismatches", "Aln Len")
        self.mtree = ttk.Treeview(t2, columns=cols, show="headings", height=20)
        for col in cols:
            w = 130 if "Seq" in col else 85
            self.mtree.heading(col, text=col)
            self.mtree.column(col, width=w, anchor="center")
        mtsb = ttk.Scrollbar(t2, orient="vertical", command=self.mtree.yview)
        self.mtree.configure(yscrollcommand=mtsb.set)
        mtsb.pack(side="right", fill="y")
        self.mtree.pack(fill="both", expand=True, padx=8, pady=8)

        # Info tab
        self.info = tk.Text(t3, wrap="word", font=FMS, bg=BG, fg=TXT,
                            relief="flat", bd=0, highlightthickness=0, state="disabled")
        ttk.Scrollbar(t3, orient="vertical", command=self.info.yview).pack(
            side="right", fill="y")
        self.info.pack(fill="both", expand=True, padx=8, pady=8)

    def _setup_tags(self):
        o = self.out
        for b, c in NUC.items():
            if b != "-":
                o.tag_configure(f"n{b}", background=c, foreground=BG, font=FMB)
        o.tag_configure("gap",   foreground="#1a4030")
        o.tag_configure("match", foreground=GLOW, font=FMB)
        o.tag_configure("mis",   foreground=RED)
        o.tag_configure("cons",  foreground=CYAN, font=FMB)
        o.tag_configure("ruler", foreground=DIM,  font=FMS)
        o.tag_configure("name",  foreground=CYAN, font=FMB)
        o.tag_configure("hdr",   foreground=GLOW,
                        font=("Courier New", 12, "bold"))
        o.tag_configure("sub",   foreground=DIM,  font=FMS)
        o.tag_configure("sm",    foreground=GLOW)
        o.tag_configure("smi",   foreground=RED)

    # ── Status ───────────────────────────────────────────────────
    def _status(self, msg, spin=False):
        self._status_lbl.config(text=msg, fg=GLOW if spin else DIM)
        if spin:
            self._pb.pack(side="left", padx=6)
            self._pb.start(8)
        else:
            self._pb.stop()
            try:
                self._pb.pack_forget()
            except Exception:
                pass

    # ── Data helpers ─────────────────────────────────────────────
    def _get_seqs(self):
        result = []
        for box in self.seq_boxes:
            raw = box.get("1.0", tk.END).strip()
            if not raw:
                continue
            if any(l.startswith(">") for l in raw.splitlines()):
                result.extend(self.aligner.parse_fasta_text(raw))
            else:
                for line in raw.splitlines():
                    if line.strip():
                        result.append((f"Seq{len(result)+1}", line.strip()))
        return result

    def _get_params(self):
        try:
            m = int(self._pvars["match"].get())
            mm = int(self._pvars["mismatch"].get())
            go = float(self._pvars["gap_open"].get())
            ge = float(self._pvars["gap_extend"].get())
            return m, mm, go, ge
        except ValueError:
            raise ValueError(
                "Scoring parameters must be numbers (gap_open/extend can be decimal).")

    # ── Run ──────────────────────────────────────────────────────
    def _run_threaded(self):
        self._status("Running…", spin=True)
        threading.Thread(target=self._run, daemon=True).start()

    def _run(self):
        try:
            seqs = self._get_seqs()
            if len(seqs) < 2:
                raise ValueError("Need at least 2 sequences.")
            m, mm, go, ge = self._get_params()
            mode = self.aln_mode.get()
            self.root.after(0, self._clear_out)

            if mode in ("global", "local"):
                self._pairwise(seqs[0], seqs[1], mode, m, mm, go, ge)
            else:
                if len(seqs) < 3:
                    raise ValueError("MSA mode needs at least 3 sequences.")
                self._msa(seqs, m, mm, go, ge)

            self.root.after(0, lambda: self._status("Done ✓"))
        except Exception as ex:
            self.root.after(0, lambda: messagebox.showerror("Error", str(ex)))
            self.root.after(0, lambda: self._status(f"Error: {ex}"))

    def _clear_out(self):
        self.out.config(state="normal")
        self.out.delete("1.0", tk.END)
        self.out.config(state="disabled")
        for i in self.mtree.get_children():
            self.mtree.delete(i)

    # ── Pairwise ─────────────────────────────────────────────────
    def _pairwise(self, t1, t2, mode, m, mm, go, ge):
        n1, s1 = t1
        n2, s2 = t2
        if mode == "global":
            a1, a2, sc = self.aligner.needleman_wunsch(s1, s2, m, mm, go, ge)
            title = "GLOBAL — NEEDLEMAN-WUNSCH  (EMBOSS Needle compatible)"
        else:
            a1, a2, sc = self.aligner.smith_waterman(s1, s2, m, mm, go, ge)
            title = "LOCAL — SMITH-WATERMAN  (EMBOSS Water compatible)"

        stats = self.aligner.compute_stats(a1, a2)
        stats["score"] = sc
        stats["seq_type"] = self.aligner.detect_sequence_type(s1)
        chunk = self.chunk.get()
        cmode = self.cmode.get()

        self.root.after(0, lambda: self._draw_pair(
            title, n1, a1, n2, a2, stats, chunk, cmode))
        self.root.after(0, lambda: self._upd_stats(stats))
        self.root.after(0, lambda: self.mtree.insert("", "end", values=(
            n1, n2, sc,
            f"{stats['identity']:.1f}%",
            f"{stats['similarity']:.1f}%",
            stats["matches"],
            stats["mismatches"],
            stats["length"],
        )))
        self.root.after(0, lambda: self._set_info(
            self._make_info([(n1, s1), (n2, s2)])))

    def _draw_pair(self, title, n1, a1, n2, a2, stats, chunk, cmode):
        o = self.out
        o.config(state="normal")

        def i(t, tg=None):
            o.insert("end", t, tg) if tg else o.insert("end", t)

        def base(ch, sym):
            if ch == "-":
                i(ch, "gap")
                return
            if cmode == "nucleotide" and ch.upper() in NUC:
                i(ch, f"n{ch.upper()}")
            else:
                i(ch, "match" if sym == "|" else "mis")

        lw = max(len(n1), len(n2)) + 2
        sym = stats["symbols"]
        i(f"── {title} ──\n", "hdr")
        i(f"Score:{stats['score']}  Identity:{stats['identity']:.2f}%  "
          f"Similarity:{stats['similarity']:.2f}%  "
          f"Matches:{stats['matches']}  Mismatches:{stats['mismatches']}  "
          f"Gaps:{stats['gaps']}  Gap%:{stats['gap_percent']:.1f}%\n\n", "sub")

        for s in range(0, len(a1), chunk):
            e = min(s + chunk, len(a1))
            c1 = a1[s:e]
            c2 = a2[s:e]
            cs = sym[s:e]
            pos1 = a1[:s].replace("-", "")
            pos2 = a2[:s].replace("-", "")
            i(f"{'':>{lw}}  {len(pos1)+1}\n", "ruler")
            i(f"{n1:<{lw}}  ", "name")
            for ch, sy in zip(c1, cs):
                base(ch, sy)
            i(f"  {len(pos1)+len(c1.replace('-',''))}\n")
            i(f"{'':>{lw}}  ")
            for sy in cs:
                i(sy, "sm") if sy == "|" else i(
                    sy, "smi") if sy == "." else i(" ")
            i("\n")
            i(f"{n2:<{lw}}  ", "name")
            for ch, sy in zip(c2, cs):
                base(ch, sy)
            i(f"  {len(pos2)+len(c2.replace('-',''))}\n\n")

        self._legend(cmode)
        o.config(state="disabled")
        self.last_out = o.get("1.0", tk.END)

    # ── MSA ──────────────────────────────────────────────────────
    def _msa(self, seqs, m, mm, go, ge):
        from itertools import combinations
        aligned, sm = self.aligner.progressive_msa(seqs)
        cons = self.aligner.msa_conservation(aligned)
        chunk = self.chunk.get()
        cmode = self.cmode.get()

        self.root.after(0, lambda: self._draw_msa(aligned, cons, chunk, cmode))

        names = [n for n, _ in seqs]
        for i, j in combinations(range(len(seqs)), 2):
            sc = sm.get((min(i, j), max(i, j)), "—")
            # Align the two raw seqs for stats
            a1t, a2t, _ = self.aligner.needleman_wunsch(
                seqs[i][1], seqs[j][1], m, mm, go, ge)
            st = self.aligner.compute_stats(a1t, a2t)
            row = (
                names[i], names[j], sc,
                f"{st['identity']:.1f}%",
                f"{st['similarity']:.1f}%",
                st["matches"],
                st["mismatches"],
                st["length"],
            )
            self.root.after(
                0, lambda r=row: self.mtree.insert("", "end", values=r))

        tl = max(len(s) for _, s in aligned)
        cv = cons.count("*")
        self.root.after(0, lambda: self._upd_stats({
            "score":       "—",
            "seq_type":    self.aligner.detect_sequence_type(seqs[0][1]),
            "length":      tl,
            "identity":    cv / tl * 100 if tl else 0,
            "similarity":  (cv + cons.count(".")) / tl * 100 if tl else 0,
            "matches":     cv,
            "mismatches":  cons.count("."),
            "gaps":        cons.count(" "),
            "gap_percent": cons.count(" ") / tl * 100 if tl else 0,
        }))
        self.root.after(0, lambda: self._set_info(self._make_info(seqs)))

    def _draw_msa(self, aligned, cons, chunk, cmode):
        o = self.out
        o.config(state="normal")

        def i(t, tg=None):
            o.insert("end", t, tg) if tg else o.insert("end", t)

        def base(ch):
            if ch == "-":
                i(ch, "gap")
                return
            i(ch, f"n{ch.upper()}") if cmode == "nucleotide" and ch.upper(
            ) in NUC else i(ch, "match")

        n = len(aligned)
        total = max(len(s) for _, s in aligned)
        lw = max(len(nm) for nm, _ in aligned) + 2
        i(f"── MULTIPLE SEQUENCE ALIGNMENT ({n} seqs, aln_len={total}) ──\n\n", "hdr")

        for s in range(0, total, chunk):
            e = min(s + chunk, total)
            i(f"{'':>{lw}}  {s+1}\n", "ruler")
            for nm, seq in aligned:
                seg = seq[s:e].ljust(e - s, "-")
                i(f"{nm:<{lw}}  ", "name")
                for ch in seg:
                    base(ch)
                i("\n")
            cs = cons[s:e].ljust(e - s)
            i(f"{'':>{lw}}  ")
            for sy in cs:
                i(sy, "cons") if sy == "*" else i(sy,
                                                  "smi") if sy == "." else i(" ")
            i("\n\n")

        self._legend(cmode)
        o.config(state="disabled")
        self.last_out = o.get("1.0", tk.END)

    def _legend(self, cmode):
        o = self.out

        def i(t, tg=None):
            o.insert("end", t, tg) if tg else o.insert("end", t)

        i("\n  Legend: ", "sub")
        if cmode == "nucleotide":
            for b in "ATGC":
                i(f" {b} ", f"n{b}")
                i("  ", "sub")
        else:
            i("█", "match")
            i(" Match  ", "sub")
            i("█", "mis")
            i(" Mismatch  ", "sub")
        i("- ", "gap")
        i("Gap  ", "sub")
        i("*", "cons")
        i(" Conserved\n", "sub")

    # ── Info tab ─────────────────────────────────────────────────
    def _make_info(self, seqs):
        lines = []
        for name, seq in seqs:
            stype = self.aligner.detect_sequence_type(seq)
            extra = ""
            if stype == "Nucleotide":
                gc = sum(1 for c in seq.upper() if c in "GC")
                extra = (f"\n  GC: {gc/len(seq)*100:.1f}%  "
                         f"A:{seq.upper().count('A')} "
                         f"T:{seq.upper().count('T')} "
                         f"G:{seq.upper().count('G')} "
                         f"C:{seq.upper().count('C')}")
            lines.append(
                f"  {name}  |  {stype}  |  {len(seq)} bp{extra}\n  {'─'*40}\n")
        return "\n".join(lines)

    def _set_info(self, text):
        self.info.config(state="normal")
        self.info.delete("1.0", tk.END)
        self.info.insert("1.0", text)
        self.info.config(state="disabled")

    # ── Utilities ─────────────────────────────────────────────────
    def _load_fasta(self):
        path = filedialog.askopenfilename(
            filetypes=[("FASTA", "*.fasta *.fa *.fna *.faa"), ("All", "*.*")])
        if not path:
            return
        try:
            seqs = self.aligner.load_fasta_file(path)
            # Auto-select mode
            if len(seqs) >= 3:
                self.aln_mode.set("msa")
            elif len(seqs) == 2:
                self.aln_mode.set("global")
            self._rebuild_input()
            # Paste into the first (only) box in the dynamic area
            self.seq_boxes[0].delete("1.0", tk.END)
            if self.aln_mode.get() in ("global", "local") and len(seqs) >= 2:
                # Two boxes available
                self._rebuild_input()
                self.seq_boxes[0].delete("1.0", tk.END)
                self.seq_boxes[0].insert(
                    "end", f">{seqs[0][0]}\n{seqs[0][1]}\n")
                self.seq_boxes[1].delete("1.0", tk.END)
                self.seq_boxes[1].insert(
                    "end", f">{seqs[1][0]}\n{seqs[1][1]}\n")
            else:
                for name, seq in seqs:
                    self.seq_boxes[0].insert("end", f">{name}\n{seq}\n")
            self._status(f"Loaded {len(seqs)} sequence(s) from file")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _clear(self):
        for box in self.seq_boxes:
            box.delete("1.0", tk.END)
        self._clear_out()
        for lbl in self._slbls.values():
            lbl.config(text="—")
        self._set_info("")
        self.last_out = ""
        self._status("Cleared")

    def _save(self):
        content = self.last_out or self.out.get("1.0", tk.END)
        if not content.strip():
            messagebox.showwarning("Empty", "Run an alignment first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".aln",
            filetypes=[("Alignment", "*.aln"), ("Text", "*.txt"), ("All", "*.*")])
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(content)
            self._status(f"Saved → {path}")


if __name__ == "__main__":
    root = tk.Tk()
    MSAApp(root)
    root.mainloop()
