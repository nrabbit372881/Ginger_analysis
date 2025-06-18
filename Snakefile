configfile: "configs/config.json"

# configure docker mounting options
# ==============================================================================
docker_mount_opt = ""
for volume in config["volumes"]:
    docker_mount_opt += "-v %s:%s:%s " % (volume["real"], volume["virtual"], volume["mode"])
    
    if volume["is_workdir"]:
        docker_mount_opt += "-w %s " % volume["virtual"]

# define query function(設計一個尋找用的函數)
# ==============================================================================
def query(**kwargs):
    for d in config["samples"]:
        for k, v in kwargs.items():
            if d[k] == v:
                return d
                
    raise ValueError("No such sample: %s" % kwargs)
    return None

# download references from NCBI dataset
# ==============================================================================
rule ncbi_dload_dehydrated: # dehydrated 模式比較容易 resume 或重跑，可rehydrate
    output:
        fullname="references/ncbi_dataset.zip"
    log:
        "logs/ncbi-dataset/ncbi_dload_dehydrated.log"
    shell:
        # line 37: dirname 是一個 bash 指令，用來從完整路徑中取出「所在的資料夾」，references/ncbi_dataset.zip的references。
        # line 41: 用我的id跑才不會超出權限
        # line 43: 本機先下載（pull）ccc/ncbi-datasets:20230926 這個 Docker image，然後根據這個 image 啟動出一個新的容器，並加上 docker_mount_opt 的 volume 掛載。
        # line 44: 我要docker執行的command
        # line 45 & 46: 1是將 stdout 寫入；2是將 stderr。  寫入兩者等同 &> {log}
        r"""
        mkdir -p $(dirname {output.fullname}) 
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name ncbi_dload_dehydrated \
            ccc/ncbi-datasets:20230926 \
                datasets download genome accession GCF_018446385.1 \
                    --include gff3,gtf,genome,seq-report \
                    --filename {output.fullname} \
                    --dehydrated \
                2> {log} \
                1> {log} 
        """

rule ncbi_rehydrate:
    input:
        fin="references/ncbi_dataset.zip"
    output:
        dout=directory("references/ncbi_dataset/data")
    log:
        "logs/ncbi-dataset/ncbi_rehydrate.log"
    shell:
        # line 62: 將 {input.fin} 這個 zip 檔解壓縮到它所在的資料夾中；如果檔案已存在則只在需要時更新（-u）；-d 後面接的是目標資料夾路徑，最後要加上要解壓的 zip 檔案路徑。
        r"""
        unzip -u -d $(dirname {input.fin}) {input.fin}

        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name ncbi_rehydrate \
            ccc/ncbi-datasets:20230926 \
                datasets rehydrate \
                    --directory $(dirname {input.fin}) \
                2> {log} \
                1> {log}
        """

rule remove_ambigous: # 不用問為甚麼
    input:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic.gtf"
    output:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic_modified.gtf"
    log:
        "logs/ncbi-dataset/remove_ambigous.log"
    shell:
        """
        grep -v F6E76_pgp044 {input} > {output}
        
        echo "Ambigous gene:" > {log}
        grep F6E76_pgp044 {input} >> {log}
        """

# QC 將目標fastq.gz檔案進行質量控制
# ==============================================================================
rule fastp:
    threads: 4
    input:
        in1=lambda wildcards: query(id=wildcards.id)["fq1"],
        in2=lambda wildcards: query(id=wildcards.id)["fq2"],
    output:
        out1="outputs/afqc/{id}/{id}_1.fastq.gz",
        out2="outputs/afqc/{id}/{id}_2.fastq.gz",
        report_json="outputs/afqc/{id}/fastp.json",
        report_html="outputs/afqc/{id}/fastp.html"
    params: #  fastp的QC參數: --detect_adapter_for_pe自動偵測paired-end的adapter；--correction對重疊區域做錯誤校正；--cut_front去除前端低品質序列；--cut_tail去除尾端低品質序列；--disable_trim_poly_g不修剪PolyG tail（適用於非 Illumina NovaSeq）
        flags=" ".join([
            "--detect_adapter_for_pe",
            "--correction",
            "--cut_front",
            "--cut_tail",
            "--disable_trim_poly_g"
        ])
    log:
        "logs/fastp/{id}.log"
    shell:
# line 123: fastp 的 image
# line 125: 加入上面那些參數
# lines 126-131: output的路徑
        """
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name fastp_{wildcards.id} \
            biocontainers/fastp:v0.20.1_cv1 \
                fastp \
                    {params.flags} \
                    --in1 {input.in1} \
                    --in2 {input.in2} \
                    --out1 {output.out1} \
                    --out2 {output.out2} \
                    --json {output.report_json} \
                    --html {output.report_html} \
                1> {log} \
                2> {log} 
        """

# The GTF file might be corrupted!
# Stop at line : NC_044775.1      RefSeq  transcript      73009   148317  .       ?       .       gene_id "F6E76_pgp044"; transcript_id "unassigned_transcript_1917"; exception "trans-splicing"; gbkey "mRNA"; gene "rps12"; locus_tag "F6E76_pgp044"; transcript_biotype "mRNA"; 
# Error Message: Strand is neither '+' nor '-'!
rule filter_gtf_strand:### 因為出現了上面錯誤，這個rule是用來過濾掉gtf檔案中strand不是+或-的行
    input:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic.gtf"
    output:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic.filtered.gtf"
    log:
        "logs/rsem/filter_gtf_strand.log"
    shell:
        """
        grep -v $'\t?\\t' {input} > {output} 2> {log}
        """

# According to the GTF file given, transcript unassigned_transcript_1917 has exons from different orientations!
### 理論上應該不能直接刪掉？
rule filter_unassigned_transcript:
    input:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic.filtered.gtf"
    output:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic.finished.gtf"
    log:
        "logs/rsem/filter_unassigned_transcript.log"
    shell:
        """
        grep -v 'unassigned_transcript_1917' {input} > {output} 2> {log}
        """


# RSEM quantification
# ==============================================================================

rule rsem_prep_ref: # 把參考基因組的fasta和gtf檔案轉成RSEM的專用index
    threads: 4
    conda:
        "envs/2023aut_ginger.yaml"
    input:
        gtf="references/ncbi_dataset/data/GCF_018446385.1/genomic.finished.gtf",
        fa="references/ncbi_dataset/data/GCF_018446385.1/GCF_018446385.1_Zo_v1.1_genomic.fna"
    output:
        dout=directory("references/rsem")
    params:
        ref_prefix="references/rsem/rsem"
    log:
        "logs/rsem/rsem_prep_ref.log"
    shell:
        # line 195: 讓 RSEM 知道 bowtie2 這個程式在哪裡
        # line 196: 指定基因註解檔的路徑
        # line 198: 指定 genome fasta 檔的路徑
        # line 199: 指定輸出的檔名前綴和路徑
        """
        mkdir -p {output.dout}

        rsem-prepare-reference \
            --num-threads {threads} \
            --bowtie2 \
            --bowtie2-path $(dirname $(which bowtie2)) \
            --gtf {input.gtf} \
            --polyA \
            {input.fa} \
            {params.ref_prefix} \
            2>&1 \
            > {log}
        """

rule rsem_cal_exp: # 把每個樣本的reads定量，算出每個基因的表現量
    threads: 12
    conda:
        "envs/2023aut_ginger.yaml"
    input:
        ref_dir="references/rsem/",
        fq1="outputs/afqc/{id}/{id}_1.fastq.gz",
        fq2="outputs/afqc/{id}/{id}_2.fastq.gz",
    output:
        dout=directory("outputs/rsem-cal-expr/{id}"),
        gene_res="outputs/rsem-cal-expr/{id}/{id}.genes.results"
    params:
        ref_prefix="references/rsem/rsem",
        output_prefix=lambda wildcards: f"outputs/rsem-cal-expr/{wildcards.id}/{wildcards.id}",
        unused_flags=" ".join(
            ["--strand-specific",]
        )
    log:
        "logs/rsem/rsem_cal_exp/{id}.log"
    shell:
        """
        mkdir -p $(dirname {output.dout})

        rsem-calculate-expression \
            --paired-end \
            --num-threads {threads} \
            --bowtie2 \
            --bowtie2-path $(dirname $(which bowtie2)) \
            --time \
            {input.fq1} \
            {input.fq2} \
            {params.ref_prefix} \
            {params.output_prefix} \
        2>&1 \
        > {log}
        """

rule rsem_dmat: # 合併所有樣本的基因表現量，把每個樣本的變成一個矩陣
    conda:
        "envs/2023aut_ginger.yaml"
    input:
        expand(
            "outputs/rsem-cal-expr/{id}/{id}.genes.results",
            id=[d["id"] for d in config["samples"]]
        )
    output:
        "outputs/rsem-dmat/all.rcounts.genes.matrix"
    log:
        "logs/rsem/rsem_dmat.log"
    shell:
        """
        rsem-generate-data-matrix {input} 1> {output} 2> {log}
        """

#先下載蛋白質需要檔案
# ==============================================================================
rule download_protein_Phytozome:
    output:
        pth_Phytozome= "references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa.gz", 
        
    log:
        Phy = "logs/protein/download_protein_by_PhytozomeV9.log", 

    shell:
        """
        mkdir -p $(dirname {output.pth_Phytozome})
        curl --cookie jgi_session=/api/sessions/cd06972c4cb4c428ce4e5fbc06d50599 --output {output.pth_Phytozome} -d '{{"ids":{{"Phytozome-167":{{"file_ids":["52b9c702166e730e43a34e56"],"top_hit":"53112a1b49607a1be0055860"}}}},"api_version":"2"}}' -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/ \

            2> {log.Phy} \
            1> {log.Phy}
        """
#之後要打開csv檔案看file_id之後再重新用curl下載一次(改file id)才是真正的txt檔案，阿記得解壓縮
rule unzip_protein_Phytozome:
    input:
        "references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa.gz"
    output:
        "references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa"
    log:
        "logs/protein/unzip_protein_by_PhytozomeV9.log"
    shell:
        """  
            unzip -o -d $(dirname {input}) {input} >> {log} 2>&1

            find $(dirname {input}) -name "*.fa.gz" -exec gunzip -c {{}} \; > {output}
            awk '/^>/{{f=1}} f' {output} > {output}.tmp
            mv {output}.tmp {output}
        """

rule download_protein_NCBI:
    output:
        pth_ncbi= "references/ginger_protein/protein.faa.zip"
    log:
        "logs/protein/download_protein_by_NCBI.log"
    shell:
        """
        mkdir -p $(dirname {output.pth_ncbi})
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name ncbi_dload_dehydrated \
            ccc/ncbi-datasets:20230926 \
                datasets download genome accession GCF_018446385.1 \
                    --include protein \
                    --filename {output.pth_ncbi} \
                    
                2> {log} \
                1> {log} 
        """
rule unzip_protein_NCBI:
    input:
        "references/ginger_protein/protein.faa.zip"
    output:
        "references/ginger_protein/ncbi_dataset/data/GCF_018446385.1/protein.faa"
    log:
        "logs/protein/unzip_protein_by_NCBI.log"
    shell:
        """
        unzip -o -d $(dirname {input}) {input}
        """

# GSEA
# ==============================================================================
# make blast database
rule makeblastdb:
    input:
        "references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa"
    output:
        "references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa.psq"
    log:
        "logs/blastp/makeblastdb.log"
    shell:
        """
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name makeblastdb \
            biocontainers/blast:v2.2.31_cv2 \
                makeblastdb \
                    -in {input} \
                    -dbtype prot \
                    2> {log} \
                    1> {log}
        """


# blast ginger proteins against arabidopsis protein
rule blastp:
    threads: 24
    input:
        query="references/ginger_protein/ncbi_dataset/data/GCF_018446385.1/protein.faa",
        db="references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa"
    output:
        "outputs/blastp/ginger_vs_ath.tsv"
    log:
        "logs/blastp/ginger_vs_ath.log"
    shell:
        """
        mkdir -p $(dirname {output})
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name blastp \
            biocontainers/blast:v2.2.31_cv2 \
                blastp \
                    -query {input.query} \
                    -db {input.db} \
                    -out {output} \
                    -outfmt 6 \
                    -num_alignments 1 \
                    -num_threads {threads} \
                    2> {log} \
                    1> {log}
        """

