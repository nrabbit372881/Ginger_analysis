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

rule remove_ambigous: # 自己找不到，不需要看是怎麼找的
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