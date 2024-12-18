#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import re
import pickle
import requests
from bs4 import BeautifulSoup
from lxml import html
import tkinter as tk
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

print("読み込み中...")
os.chdir(os.path.dirname(os.path.abspath(__file__)))

plt.rcParams['font.family'] = 'Arial'  # ここでArialフォントを指定
plt.rcParams['font.size']           = 11    #　フォントサイズ
plt.rcParams['lines.linewidth']     = 1.5   # プロットの線の太さ
plt.rcParams['axes.linewidth']      = 1.5   # 枠線（スパイン）の太さ
plt.rcParams['xtick.major.width']   = 1.5   # x軸目盛りの線の太さ
plt.rcParams['ytick.major.width']   = 1.5   # y軸目盛りの線の太さ
gene_AGI=None

genes_data = pd.read_csv("Araport11_GTF_genes_transposons.csv") # GTFcsvの読み込み (遺伝子領域情報)
fasta_file = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"  # DNA配列の読み込み
gene_aliases = pd.read_csv('gene_aliases_20230630.csv')

dna_sequences ={}
for record in SeqIO.parse(fasta_file, "fasta"):
    dna_sequences[str(record.id)] = record.seq


class GeneCoveragePlotter:
    def __init__(self, genes_data):
        self.genes_data = genes_data
    

    def convert_genesymbol(self, AGI, filename = 'AGItoSymbol.pkl'):
        existing_dic ={}
        with open(filename, 'rb') as file:
            existing_dic = pickle.load(file)
        if AGI not in existing_dic.keys():
            response = requests.get(f"https://bar.utoronto.ca/thalemine/portal.do?externalids={AGI}")
            soup = BeautifulSoup(response.content, "html.parser")
            lxml_data = html.fromstring(str(soup))
            symbol  = lxml_data.xpath('//*[@id="object_header"]/a/h1/strong')[0].text   # XPathが指定できる
            if symbol == AGI:   #symbolなかった場合
                try: symbol = list(gene_aliases[gene_aliases["locus_name"] == AGI]["symbol"])[0]    #gene_aliasesからも検索する
                except: symbol == None
            
            existing_dic[AGI] = symbol
            with open(filename, 'wb') as file:
                pickle.dump(existing_dic, file)
                print("updated pkl")
            
        else:
            symbol = existing_dic[AGI]
        return symbol


    def search_nearby_region(self, start, chrom, end, gene_AGI, upper, lower):
        """
        AGIが優先的に検索される、hitしないor AGI=Noneのとき、satr, end positionで検索
        遺伝子にgene領域がない場合（non coding region）は遺伝子名を表示しない
        """
        if gene_AGI is not None:
            strand  = self.genes_data[self.genes_data['gene_id'].str.contains(gene_AGI)]["strand"].unique()[0]
            chrom   = self.genes_data[self.genes_data['gene_id'].str.contains(gene_AGI)]["chromosome"].unique()[0][-1]
            start   = self.genes_data[self.genes_data['gene_id'].str.contains(gene_AGI)]["feature_start"].min()
            end     = self.genes_data[self.genes_data['gene_id'].str.contains(gene_AGI)]["feature_end"].max()
            if strand == '-': return str(chrom) ,int(start- lower), int(end + upper)
            else: return str(chrom), int(start- upper), int(end + lower)
        else: return str(chrom), int(start), int(end)


    def plot_coverage(self, ax, bam_file, chrom, start, end, color, top):
        # BAMファイルを開く
        samfile = pysam.AlignmentFile(bam_file, "rb")
        # Calculate coverage
        coverage = [0] * (end - start)
        for read in samfile.fetch(chrom, start, end):
            if not read.is_unmapped:
                for i in range(max(read.reference_start, start), min(read.reference_end, end)):
                    coverage[i - start] += 1

        # Plot coverage
        window_size = 2  # Moving average window size
        moving_avg = np.convolve(coverage, np.ones(window_size) / window_size, mode='same')
        ax.plot(range(start, end), moving_avg, color=color, linewidth=0.5, label="Moving Average")
        ax.fill_between(range(start, end), moving_avg, color=color, alpha=0.8, label="Coverage")

        title = bam_file.split("/")[-1].split(".")[0]
        ax.set_ylabel(f"{title}", rotation=90, labelpad=5, fontsize =8)  # rotation=0で横向き、labelpadで軸からの距離を調整

        ax.set_ylim(bottom=0)
        ax.set_ylim(top=top)
        ax.set_xlim(start, end)
        ax.plot([end-1000, end], [top*0.95, top*0.95], color="black", linestyle="-", lw=0.75)
        # ax.set_ylabel("Coverage")
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=True, useMathText=True))
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        # ax.yaxis.set_major_locator(ticker.MultipleLocator(top/2))
        samfile.close()


    def plot_gene_region(self, ax, feature_type, feature_start, feature_end, strand, gene_id, gene_AGI):
        if gene_AGI == gene_id: 
            color="black"
            zorder=2
            fontsize_adj = 0
            symbol = self.convert_genesymbol(gene_AGI)
            gene_id = gene_id + " (" + symbol + ")"
            
        else: 
            color = "gray"
            zorder=1
            fontsize_adj = 1.5
        
        xlim_start, xlim_end = ax.get_xlim()
        center_point = (feature_start + feature_end) / 2

        if feature_type == "gene":
            ax.plot([feature_start, feature_end], [0, 0], color=color, linestyle="-", zorder=zorder)    #geneの線
            if center_point > xlim_start and center_point < xlim_end:
                ax.text(center_point, 2, gene_id, fontsize=10 - fontsize_adj, 
                        ha='center', va='bottom', color=color, zorder=zorder) #AGIコードのplot
        else:
            rect = patches.Rectangle((feature_start, -1),  feature_end-feature_start, 2, color=color)   #exonのbox
            ax.add_patch(rect)
            if (feature_end - feature_start) >= 100:
                mid_point = (feature_start + feature_end) / 2
                if xlim_start < center_point < xlim_end:
                    arrow = '>' if strand == '+' else '<'
                    ax.text(mid_point, -1, arrow, fontsize=10, ha='center',  color='white', fontweight='bold', zorder=zorder)


    def plot_genebody(self, ax, chrom, start, end, gene_AGI):   #plot gene body
        # ax.set_xlabel("Genomic Position (bp)")
        ax.set_xlim(start, end)
        ax.set_ylim(-3,5)
        ax.set_yticks([])  # Hide Y-axis ticks

        genes_data_nearby = self.genes_data[
            (self.genes_data["chromosome"] == f"Chr{chrom}") &
            (self.genes_data["feature_start"] <= end) &
            (self.genes_data["feature_end"] >= start)
        ].copy()

        # CDSを持たないタイプのgeneの場合は可視化しない
        if "CDS" in genes_data_nearby["feature_type"].unique():
            vis_target = "CDS"
        else: vis_target = "X"  #"exon"
            

        for _, gene_body in genes_data_nearby.iterrows():
            feature_type = gene_body["feature_type"]
            feature_start = gene_body["feature_start"]
            feature_end = gene_body["feature_end"]
            strand = gene_body["strand"]
            gene_id = gene_body["gene_id"]

            if feature_type =="gene" or feature_type == vis_target:
                self.plot_gene_region(ax, feature_type, feature_start, feature_end, strand, gene_id, gene_AGI)

        # Hide ax spines
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
        ax.tick_params(bottom=False, left=False, right=False, top=False)
        ax.set_xlim(start, end)
    

    def find_motif(self, chrom, start, end, motif_list):
        motifs = {}
        sequence = str(dna_sequences[chrom][start:end])
        
        for motif in motif_list:
            forward = str(Seq(motif))
            reverse = str(Seq(motif).reverse_complement())
            # 順方向検索
            forward_positions = [m.span() for m in re.finditer(forward, sequence)]
            if forward_positions:
                motifs[forward] = forward_positions
            else: motifs[forward] = []
            # 逆方向検索
            reverse_positions = [m.span() for m in re.finditer(reverse, sequence)]
            if reverse_positions:
                motifs[reverse] = reverse_positions
            else: motifs[reverse] = []
        return motifs


    def plot_motif(self, ax, start, end, motifs):
        ax.set_xlim(start, end)
        ax.set_ylim(0, len(motifs))
        ax.set_yticks([len(motifs) - i - 0.5 for i in range(len(motifs))])  # y位置を動的に設定
        ax.set_yticklabels(list(motifs.keys()), fontsize=6)

        for i, (label, intervals) in enumerate(motifs.items()):
            y_center = len(motifs) - i - 0.5  # 各カテゴリのy中心位置
            for (x_start, x_end) in intervals:
                rect = patches.Rectangle((start + x_start, y_center - 0.25),  x_end - x_start, 0.5, color='black')
                ax.add_patch(rect)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
        ax.tick_params(bottom=False, left=False, right=False, top=False)
        

    def make_plot_area(self, num_axes, fig_height, fig_width, fixed_height=0.6):
        """
        サブプロットを作成し、1つは高さを固定し、残りは fig の大きさに合わせて可変にする

        Args:
            num_axes (int): サブプロットの総数
            fixed_height (float): 固定するサブプロットの高さ
            fig_height (float): fig の全体の高さ (インチ)

        Returns:
            tuple: (fig, axes)
        """
        variable_height = (fig_height - 2*fixed_height) / (num_axes - 2) if num_axes > 2 else 0 # 可変部分の高さを計算
        height_ratios = [variable_height] * (num_axes - 2) + [fixed_height] + [fixed_height] # 高さの比率リストを作成
        
        gs = gridspec.GridSpec(num_axes, 1, height_ratios=height_ratios)
        fig = plt.figure(figsize=(fig_width, fig_height))
        axes = []
        for i in range(num_axes):
            axes.append(fig.add_subplot(gs[i]))
        return fig, axes
    

    def plot_multiple_coverages(self, bam_files, chrom, start, end, gene_AGI, upper, lower, top, color, fig_height, fig_width, out_path, motif_list, find_motif, save):
        new_chrom, new_start, new_end = self.search_nearby_region(chrom, start, end, gene_AGI, upper, lower)
        
        #figのサイズ比率
        num_axes = len(bam_files) + 2   #bam files + GeneBody + MotifPlot
        fig, axes = self.make_plot_area(num_axes, fig_height, fig_width)      

        # 異なるBAMファイルのカバレッジをプロット
        for i, bam_file in enumerate(bam_files):
            self.plot_coverage(axes[i], bam_file, new_chrom, new_start, new_end, color, top)
            bbox = axes[i].get_position()
            axes[i].set_position([bbox.x0, bbox.y0, bbox.width, bbox.height * 1.05])

        self.plot_genebody(axes[i+1], new_chrom, new_start, new_end, gene_AGI)
        bbox = axes[i+1].get_position()
        if find_motif:
            motifs = self.find_motif(new_chrom, new_start, new_end, motif_list)
            self.plot_motif(axes[-1], new_start, new_end, motifs) 
        else:
            for spine in axes[-1].spines.values():
                spine.set_visible(False)
            axes[-1].tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
            axes[-1].tick_params(bottom=False, left=False, right=False, top=False)
        # axes[i+1].set_position([bbox.x0, bbox.y0 * 0.8, bbox.width, bbox.height])
        if save: 
            if gene_AGI is not None: plt.savefig(f"{out_path}/{gene_AGI}_peak.png", dpi=800, transparent=True)
            else: plt.savefig(f"{out_path}/chr{new_chrom}_{new_start}_{new_end}_peak.png", dpi=800, transparent=True)
        return fig
        # plt.show()


# %%
# App GUI
class GUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Gene Coverage Plotter")

        # 左側の入力フィールド
        self.input_frame = tk.Frame(root)
        self.input_frame.pack(side=tk.LEFT, padx=10, pady=10)

        self.create_input_fields()

        # 右側のプロット表示エリア
        self.plot_frame = tk.Frame(root)
        self.plot_frame.pack(side=tk.RIGHT, padx=10, pady=10)

        self.canvas = None

    def create_input_fields(self):
        self.labels = {
            "bam_files": "BAM Files",
            "chrom": "Chromosome",
            "start": "Start Position",
            "end": "End Position",
            "gene_AGI": "Gene AGI",
            "upper": "GenomeRange Upper bp",
            "lower": "GenomeRange Lower bp",
            "top": "Y-axis Max Value",
            "color": "Color",
            "fig_height": "Figure Height",
            "fig_width": "Figure Width",
            "motif_list": "Find Motif",
        }

        self.default_values = {
            "bam_files": "",  # ファイルは選択時に設定
            "chrom": "1",
            "start": "1000",
            "end": "2000",
            "gene_AGI": "AT1G01010",
            "upper": "2000",
            "lower": "100",
            "top": "120",
            "color": "blue",
            "fig_height": "6",
            "fig_width": "10",
            "motif_list": "",
        }

        self.entries = {}

        # inputファイル選択ボタン
        self.file_button = tk.Button(
            self.input_frame, text="Select BAM Files", command=self.select_files
        )
        self.file_button.grid(row=0, column=0, columnspan=2)

        for idx, (key, label) in enumerate(self.labels.items()):
            tk.Label(self.input_frame, text=label).grid(row=idx+1, column=0, sticky=tk.W)
            entry = tk.Entry(self.input_frame)
            entry.insert(0, self.default_values[key])  # デフォルト値を設定
            entry.grid(row=idx+1, column=1)
            self.entries[key] = entry
        
        # チェックボックス (save)
        self.save_var = tk.BooleanVar(value=False)
        self.save_check = tk.Checkbutton(self.input_frame, text="Save plot", variable=self.save_var)
        self.save_check.grid(row=len(self.labels)+2, column=0, columnspan=2, sticky=tk.W)

        # チェックボックス (plot_motif)
        self.motif_var = tk.BooleanVar(value=False)
        self.motif_check = tk.Checkbutton(self.input_frame, text="Find Motif", variable=self.motif_var)
        self.motif_check.grid(row=len(self.labels)+3, column=0, columnspan=2, sticky=tk.W)

        # outputフォルダ選択ボタン
        self.output_button = tk.Button(
            self.input_frame, text="Select output folder", command=self.output_dir
        )
        self.output_button.grid(row=len(self.labels)+4, column=0, columnspan=2)

        tk.Label(self.input_frame, text="Output Path").grid(row=len(self.labels)+5, column=0, sticky=tk.W)
        entry = tk.Entry(self.input_frame)
        entry.insert(0, "")  # デフォルト値を設定
        entry.grid(row=len(self.labels)+5, column=1)
        self.entries["output_path"] = entry

        # プロットボタン
        self.plot_button = tk.Button(
            self.input_frame, text="Plot", command=self.plot_graph
        )
        self.plot_button.grid(row=len(self.labels) +6, column=0, columnspan=2)

    def select_files(self):
        file_paths = filedialog.askopenfilenames(
            title="Select BAM Files",
            filetypes=[("All Files", "*.*")]
        )
        if file_paths:
            self.entries["bam_files"].delete(0, tk.END)
            self.entries["bam_files"].insert(0, ",".join(file_paths))  # 複数ファイルをカンマで区切って設定

    def output_dir(self):
        iDir = os.path.abspath(os.path.dirname(__file__))
        iDirPath = filedialog.askdirectory(initialdir = iDir) 

        if iDirPath:
            self.entries["output_path"].delete(0, tk.END)
            self.entries["output_path"].insert(0, "".join(iDirPath))  # 複数ファイルをカンマで区切って設定

    def plot_graph(self):
        # 入力値を取得
        try:
            bam_files = self.entries["bam_files"].get().split(",")  # カンマで分割
            chrom = self.entries["chrom"].get()
            start = int(self.entries["start"].get())
            end = int(self.entries["end"].get())
            gene_AGI = self.entries["gene_AGI"].get()
            upper = float(self.entries["upper"].get())
            lower = float(self.entries["lower"].get())
            top = float(self.entries["top"].get())
            color = self.entries["color"].get()
            fig_height = float(self.entries["fig_height"].get())
            fig_width = float(self.entries["fig_width"].get())
            motif_list = [item for item in ",".join(self.entries["motif_list"].get().split(" ")).split(",") if len(item) > 0]
            if self.entries["output_path"].get() is not None:
                output_path = self.entries["output_path"].get()
            else: output_path = None
            save = self.save_var.get()
            find_motif = self.motif_var.get()

            plotter = GeneCoveragePlotter(genes_data)
            fig = plotter.plot_multiple_coverages(
                bam_files, chrom, start, end, gene_AGI, upper, lower, top, 
                color, fig_height, fig_width, output_path, motif_list, find_motif=find_motif, save=save,
            )

            self.update_plot(fig)

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_plot(self, fig):
        # 既存のキャンバスを削除
        if self.canvas:
            self.canvas.get_tk_widget().destroy()

        # 新しいFigureを表示
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()


if __name__ == "__main__":
    root = tk.Tk()
    app = GUI(root)
    root.mainloop()