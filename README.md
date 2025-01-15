# ChIP-seq Peak Visualize アプリの実行方法

## 必要要件

- Python 3.10 推奨
- 必要な Python ライブラリは `requirements.txt` に記載されています。

## 実行手順（ターミナル）

1. **実行環境を構築する**

   ```bash
   git clone https://github.com/PM2951/ChIP-seq_PeakVis.git
   cd ChIP-seq_PeakVis
   pip install -r requirements.txt
   curl https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -o Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
   gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
   unzip Araport11_GTF_genes_transposons.csv.zip
   ```

   エラーが出る場合は仮想環境を構築後に実行してください。

2. **アプリを実行する**

   以下のコマンドでアプリケーションを起動します。

   ```bash
   python ChIP-seq_Peak_visualize_app.py
   ```

   もしくは、
   
   ```bash
   python ChIP-seq_PeakVis/ChIP-seq_Peak_visualize_app.py
   ```
   

4. **アプリケーションの終了**

   アプリの終了ボタン
   または、ターミナル上で cntrol + Cボタン を押し、終了させる。


## 🔧 パラメーターについて
### **入力ファイル**
- **説明:** 可視化したい`.sort.bam`ファイルを1つ以上選択します。
- **注意:** 対応する`.sort.bam.bai`インデックスファイルが同じディレクトリに存在する必要があります。

### **モチーフ配列検索**
- **説明:** `.sort.bam`ファイル内で検索したいモチーフ配列を入力します。
- **フォーマット:** モチーフ配列をカンマ（`,`）区切りで入力します。例: `ATG, TATA, CCGG`。
- チェックボックスをクリックすることでplotされます。

## 📊 可視化のカスタマイズ
アプリケーション内で以下のパラメーターを調整できます。
- **Gene AGI:** Arabidopsis Thaliana のAGI codeを大文字で入力してください。例) AT1G01010
- **GenomeRange Upper/Lower bp:** 遺伝子領域からの距離を設定します。strandを軸に上流をUpper, 下流をLowerとしています。
- **プロットの幅と高さ:** 生成されるプロットのサイズを設定します。
- **Color:** カバレッジプロットの色をカスタマイズします。black, #ffffff などPythonで用いられるカラーコードが対応しています。
- **Y-axis Max Value:** プロットのy軸の範囲を設定します。
- **Prefix:** 出力ファイルの名前に付加するプレフィックスを指定します。

## 📂 出力
可視化されたプロットはPNG画像（800 dpi, 背景透過）として保存されます。出力パスとファイル名のプレフィックスはアプリケーション内で指定できます。

## 📄 ライセンス
このプロジェクトは名古屋大学理学研究科多田研究室によって作成されました。詳細についてはご連絡ください。

https://www.gene.nagoya-u.ac.jp/

