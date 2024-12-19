# ChIP-seq Peak Visualize アプリの実行方法

## 必要要件

- Python 3.7 以上
- 必要な Python ライブラリは `requirements.txt` に記載されています。

## 実行手順（ターミナル）

1. **リポジトリをクローンする**

   リポジトリをローカル環境にクローンします。

   ```bash
   git clone https://github.com/PM2951/ChIP-seq_PeakVis.git
   cd ChIP-seq_PeakVis
   ```

   もしくは、zipファイルをダウンロードして解凍してください。

2. **必要なライブラリをインストールする**

   必要な Python ライブラリをインストールします。

   ```bash
   pip install -r requirements.txt
   ```

   Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gzをダウンロードします。（Windows の場合: curl　ではなく wget　 が実行可能な場合もあります。)

   ```bash
   curl https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -o Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
   ```

   圧縮データを解凍します。

   ```bash
   gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
   unzip Araport11_GTF_genes_transposons.csv.zip
   ```

3. **アプリを実行する**

   以下のコマンドでアプリケーションを起動します。

   ```bash
   python ChIP-seq_Peak_visualize_app.py
   ```

   もしくは、
   
   ```bash
   python ChIP-seq_PeakVis/ChIP-seq_Peak_visualize_app.py
   ```
   

4. **アプリケーションの終了**

    control + qボタン を押し、終了させる。


### パラメーター
   可視化したい .sort.bamファイルは一度に複数選択できます。
   ＊ .sort.bam.baiが同階層に存在する必要があります。

   検索したいMotif配列は カンマ区切り(,)で入力してください。

